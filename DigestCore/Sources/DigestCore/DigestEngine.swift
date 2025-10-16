import Foundation

public struct DigestOptions {
    public let circular: Bool
    public let returnSequences: Bool
    public init(circular: Bool = false, returnSequences: Bool = false) {
        self.circular = circular
        self.returnSequences = returnSequences
    }
}

public final class DigestEngine {
    let sequence: String
    let enzymes: [Enzyme]
    // Performance optimization: cached sequence as character array
    let sequenceChars: [Character]
    
    public init(sequence: String, enzymes: [Enzyme]) {
        self.sequence = sequence.uppercased()
        self.enzymes = enzymes
        self.sequenceChars = Array(self.sequence)
    }

    public func cutSites() -> [CutSite] {
        var sites: [CutSite] = []
        for e in enzymes {
            let hits = findAllMatches(sequenceChars: sequenceChars, motifChars: e.siteChars)
            for start in hits {
                // Use only cutIndexTop for the cut position
                // cutIndexBottom is used only for overhang calculation
                if let top = e.cutIndexTop {
                    sites.append(.init(enzyme: e, position: start + top, strand: .top))
                }
            }
        }
        let uniq = Array(Set(sites))
        return uniq.sorted { $0.position < $1.position }
    }

    public func digest(options: DigestOptions) -> [Fragment] {
        let n = sequence.count
        let allCutSites = cutSites()
        let cuts = allCutSites.map(\.position).sorted()
        
        guard !cuts.isEmpty else {
            // No cuts: return entire sequence with natural termini
            let naturalEnd = EndInfo(overhangType: .blunt, overhangSeq5to3: nil, sourceEnzyme: nil, overhangLen: 0)
            return [Fragment(
                start: 0, end: n, length: n,
                leftEnd: naturalEnd,
                rightEnd: naturalEnd,
                sequence: options.returnSequences ? sequence : nil
            )]
        }

        var frags: [Fragment] = []
        if options.circular {
            let extended = cuts + [cuts.first! + n]
            // Create position → CutSite lookup (pick first if multiple at same position)
            let siteByPos = Dictionary(allCutSites.map { ($0.position, $0) }, uniquingKeysWith: { first, _ in first })
            
            for k in 0..<extended.count-1 {
                let a = extended[k]
                let b = extended[k+1]
                let start = a % n
                let end = b % n
                let len = end >= start ? (end - start) : (n - start + end)

                let left = endInfoForBoundary(boundaryPos: start,
                                               cutAtPos: siteByPos[a % n],
                                               circular: true)
                let right = endInfoForBoundary(boundaryPos: end,
                                                cutAtPos: siteByPos[b % n],
                                                circular: true)

                frags.append(Fragment(
                    start: start, end: end, length: len,
                    leftEnd: left, rightEnd: right,
                    sequence: options.returnSequences ? sliceCircular(sequenceChars, start, len) : nil
                ))
            }
        } else {
            let indices = [0] + cuts + [n]
            let siteByPos = Dictionary(allCutSites.map { ($0.position, $0) }, uniquingKeysWith: { first, _ in first })
            
            for i in 0..<indices.count-1 {
                let start = indices[i]
                let end = indices[i+1]
                let leftCut = (i == 0) ? nil : siteByPos[start]
                let rightCut = (i == indices.count-2) ? nil : siteByPos[end]
                
                let left = endInfoForBoundary(boundaryPos: start, cutAtPos: leftCut, circular: false)
                let right = endInfoForBoundary(boundaryPos: end, cutAtPos: rightCut, circular: false)

                frags.append(Fragment(
                    start: start, end: end, length: end - start,
                    leftEnd: left, rightEnd: right,
                    sequence: options.returnSequences ? String(sequenceChars[start..<end]) : nil
                ))
            }
        }
        return frags
    }
    
    private func endInfoForBoundary(boundaryPos: Int,
                                    cutAtPos: CutSite?,
                                    circular: Bool) -> EndInfo {
        guard let cs = cutAtPos else {
            // Natural terminus (linear sequence end)
            return EndInfo(overhangType: .blunt, overhangSeq5to3: nil, sourceEnzyme: nil, overhangLen: 0)
        }
        
        let e = cs.enzyme
        let typ = e.overhangType ?? .unknown
        
        // Fall back: infer type if missing
        let overhangType: OverhangType = (typ == .unknown)
            ? inferType(siteLen: e.site.count, cutIndexTop: e.cutIndexTop)
            : typ

        let k = (e.cutIndexTop == nil) ? 0 : overhangLength(site: e.site, cutIndexTop: e.cutIndexTop!)
        let sticky = (e.cutIndexTop == nil) ? "" :
            canonicalSticky(seq: sequence,
                            cutPos: boundaryPos,
                            recognitionSite: e.site,
                            cutIndexTop: e.cutIndexTop!,
                            overhangType: overhangType,
                            circular: circular)
        
        return EndInfo(overhangType: overhangType,
                       overhangSeq5to3: sticky.isEmpty ? nil : sticky,
                       sourceEnzyme: e.name,
                       overhangLen: k)
    }
    
    private func inferType(siteLen: Int, cutIndexTop: Int?) -> OverhangType {
        guard let c = cutIndexTop else { return .blunt }
        let d = 2*c - siteLen
        if d == 0 { return .blunt }
        // Sign convention: positive → 5' overhang (matches Python's examples like EcoRI)
        return d > 0 ? .fivePrime : .threePrime
    }
}

private func sliceCircular(_ arr: [Character], _ start: Int, _ length: Int) -> String {
    let n = arr.count
    return String((0..<length).map { arr[(start + $0) % n] })
}


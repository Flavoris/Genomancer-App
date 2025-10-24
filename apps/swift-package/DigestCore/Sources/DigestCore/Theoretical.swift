import Foundation

public func theoreticalTemplate(site: String, cutIndexTop: Int, typ: OverhangType, k: Int) -> String {
    guard k > 0 else { return "" }
    let s = site.uppercased()
    if typ == .fivePrime {
        // bases AFTER the cut within/near the site
        if cutIndexTop + k <= s.count {
            let arr = Array(s); return String(arr[cutIndexTop..<cutIndexTop+k])
        } else {
            // Type IIS: extend with "N" beyond site bounds
            let arr = Array(s)
            let inside = max(0, s.count - cutIndexTop)
            let head = String(arr[cutIndexTop..<s.count])
            return head + String(repeating: "N", count: k - inside)
        }
    } else if typ == .threePrime {
        if cutIndexTop - k >= 0 {
            let arr = Array(s); return String(arr[cutIndexTop-k..<cutIndexTop])
        } else {
            let arr = Array(s)
            let inside = max(0, cutIndexTop)
            let tail = String(arr[0..<cutIndexTop])
            return String(repeating: "N", count: k - inside) + tail
        }
    } else { return "" }
}

public struct TheoreticalEnd {
    public let enzyme: String
    public let type: OverhangType
    public let k: Int
    public let template: String
    public let palindromic: Bool
}

public func theoreticalEnd(enzyme e: Enzyme) -> TheoreticalEnd {
    let t = e.overhangType ?? .unknown
    if t == .blunt { return .init(enzyme: e.name, type: .blunt, k: 0, template: "", palindromic: false) }
    guard let ci = e.cutIndexTop else { return .init(enzyme: e.name, type: .blunt, k: 0, template: "", palindromic: false) }
    let k = overhangLength(site: e.site, cutIndexTop: ci)
    let tpl = theoreticalTemplate(site: e.site, cutIndexTop: ci, typ: t, k: k).uppercased()
    let pal = (k > 0) && (tpl == revcomp(tpl))
    return .init(enzyme: e.name, type: t, k: k, template: tpl, palindromic: pal)
}

public func compatibleTheoretical(a: TheoreticalEnd, b: TheoreticalEnd,
                                  includeBlunt: Bool = false,
                                  minOverhang: Int = 1) -> Bool {
    if a.k == 0 && b.k == 0 { return includeBlunt }
    if a.k == 0 || b.k == 0 { return false }
    if a.type != b.type { return false }
    if a.k < minOverhang || b.k < minOverhang { return false }
    if a.k != b.k { return false }
    
    // For sticky ends to be compatible, template A should match revcomp(template B)
    // This matches the logic in compatibleConcrete: sa == revcomp(sb)
    let templateA = a.template.uppercased()
    let templateB = b.template.uppercased()
    let templateBRevComp = revcomp(templateB).uppercased()
    
    // Check if templateA matches templateBRevComp using IUPAC character matching
    let arrA = Array(templateA)
    let arrB = Array(templateBRevComp)
    guard arrA.count == arrB.count else { return false }
    
    for i in 0..<arrA.count {
        guard let setA = IUPAC.map[arrA[i]], let setB = IUPAC.map[arrB[i]] else { return false }
        // Check if the sets have any overlap (can represent the same base)
        if setA.intersection(setB).isEmpty { return false }
    }
    return true
}


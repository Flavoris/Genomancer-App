import SwiftUI
import DigestCore

struct CutSiteDetailView: View {
    let cutSite: CutSite
    let index: Int
    let sequence: String
    
    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 20) {
                // Cut Site Summary Card
                cutSiteSummaryCard
                
                // Cut Site Ends Information
                cutSiteEndsCard
            }
            .padding()
        }
        .background(Color.genomancerBackground)
        .navigationTitle("Cut Site \(index + 1)")
        .navigationBarTitleDisplayMode(.inline)
        .toolbarColorScheme(.dark, for: .navigationBar)
    }
    
    private var cutSiteSummaryCard: some View {
        VStack(alignment: .leading, spacing: 12) {
            Text("Cut Site Information")
                .font(.headline)
                .foregroundColor(.primary)
            
            Divider()
            
            HStack {
                Text("Enzyme:")
                    .foregroundColor(.secondary)
                Spacer()
                Text(cutSite.enzyme.name)
                    .font(.headline)
                    .foregroundColor(.genomancerAccent)
            }
            
            HStack {
                Text("Position:")
                    .foregroundColor(.secondary)
                Spacer()
                Text("\(cutSite.position) bp")
                    .font(.system(.body, design: .monospaced))
                    .foregroundColor(.primary)
            }
            
            HStack {
                Text("Recognition Site:")
                    .foregroundColor(.secondary)
                Spacer()
                Text(cutSite.enzyme.site)
                    .font(.system(.body, design: .monospaced))
                    .foregroundColor(.genomancerInfo)
            }
        }
        .padding()
        .background(Color(.systemBackground))
        .cornerRadius(12)
    }
    
    private var cutSiteEndsCard: some View {
        VStack(alignment: .leading, spacing: 16) {
            Text("Cut Site Ends")
                .font(.headline)
                .foregroundColor(.primary)
            
            Divider()
            
            // Left side of cut
            VStack(alignment: .leading, spacing: 8) {
                HStack {
                    Text("Left Side")
                        .font(.subheadline)
                        .fontWeight(.semibold)
                        .foregroundColor(.genomancerAccent)
                    Spacer()
                }
                
                cutEndDetails(isLeftEnd: true)
            }
            .padding()
            .background(Color(.systemGray6))
            .cornerRadius(8)
            
            // Right side of cut
            VStack(alignment: .leading, spacing: 8) {
                HStack {
                    Text("Right Side")
                        .font(.subheadline)
                        .fontWeight(.semibold)
                        .foregroundColor(.genomancerAccent)
                    Spacer()
                }
                
                cutEndDetails(isLeftEnd: false)
            }
            .padding()
            .background(Color(.systemGray6))
            .cornerRadius(8)
        }
        .padding()
        .background(Color(.systemBackground))
        .cornerRadius(12)
    }
    
    private func cutEndDetails(isLeftEnd: Bool) -> some View {
        VStack(alignment: .leading, spacing: 8) {
            HStack {
                Text("Enzyme:")
                    .foregroundColor(.secondary)
                    .font(.caption)
                Text(cutSite.enzyme.name)
                    .font(.caption)
                    .fontWeight(.medium)
                    .foregroundColor(.primary)
            }
            
            HStack {
                Text("Type:")
                    .foregroundColor(.secondary)
                    .font(.caption)
                Text(overhangDescription)
                    .font(.caption)
                    .foregroundColor(.primary)
            }
            
            if let visualization = createCutSiteVisualization(isLeftEnd: isLeftEnd) {
                VStack(alignment: .leading, spacing: 4) {
                    Text("Overhang:")
                        .foregroundColor(.secondary)
                        .font(.caption)
                    
                    VStack(alignment: .leading, spacing: 2) {
                        Text(visualization.topStrand)
                            .font(.system(.caption, design: .monospaced))
                            .foregroundColor(.primary)
                        Text(visualization.bottomStrand)
                            .font(.system(.caption, design: .monospaced))
                            .foregroundColor(.primary)
                    }
                    .padding(8)
                    .background(Color.genomancerAccent.opacity(0.1))
                    .cornerRadius(6)
                }
            }
        }
    }
    
    private var overhangDescription: String {
        if let overhangType = cutSite.enzyme.overhangType {
            switch overhangType {
            case .blunt:
                return "Blunt"
            case .fivePrime:
                return "5' overhang"
            case .threePrime:
                return "3' overhang"
            case .unknown:
                return "Unknown"
            }
        }
        
        // Fallback: calculate from cut indices
        guard let cutTop = cutSite.enzyme.cutIndexTop,
              let cutBottom = cutSite.enzyme.cutIndexBottom else {
            return "Unknown"
        }
        
        if cutTop == cutBottom {
            return "Blunt"
        } else if cutTop < cutBottom {
            return "5' overhang"
        } else {
            return "3' overhang"
        }
    }
    
    // MARK: - DNA Visualization
    
    private struct DNAVisualization {
        let topStrand: String
        let bottomStrand: String
    }
    
    private func complement(_ base: Character) -> Character {
        switch base.uppercased().first {
        case "A": return "T"
        case "T": return "A"
        case "G": return "C"
        case "C": return "G"
        default: return "N"
        }
    }
    
    private func complementSequence(_ seq: String) -> String {
        return String(seq.map { complement($0) })
    }
    
    private func createCutSiteVisualization(isLeftEnd: Bool) -> DNAVisualization? {
        let site = cutSite.enzyme.site
        guard let cutTop = cutSite.enzyme.cutIndexTop,
              let cutBottom = cutSite.enzyme.cutIndexBottom else {
            return nil
        }
        
        let position = cutSite.position  // This is the cut position (0-indexed)
        let siteLength = site.count
        
        // Calculate where the recognition site starts (going backwards from cut position)
        let siteStart = position - cutTop
        
        // Validate bounds
        guard siteStart >= 0 && siteStart + siteLength <= sequence.count else {
            return nil
        }
        
        // Get context bases for display
        let contextBases = 3
        
        // Determine overhang size and type
        let overhangSize = abs(cutTop - cutBottom)
        
        var topStrand: String
        var bottomStrand: String
        
        if cutTop == cutBottom {
            // Blunt end - cut at same position on both strands
            let cutPos = position  // This is already the cut position
            
            if isLeftEnd {
                // Left side: show bases before the cut
                let start = max(0, cutPos - contextBases)
                let end = cutPos
                let leftSeq = extractSequence(start: start, end: end)
                let leftComp = complementSequence(leftSeq)
                
                topStrand = "\(leftSeq) -5'"
                bottomStrand = "\(leftComp) -3'"
            } else {
                // Right side: show bases after the cut
                let start = cutPos
                let end = min(sequence.count, cutPos + contextBases)
                let rightSeq = extractSequence(start: start, end: end)
                let rightComp = complementSequence(rightSeq)
                
                topStrand = "\(rightSeq) -5'"
                bottomStrand = "\(rightComp) -3'"
            }
        } else if cutTop < cutBottom {
            // 5' overhang - top strand cuts earlier, bottom strand cuts later
            
            if isLeftEnd {
                // Left side: top strand is shorter, bottom strand extends with overhang
                // The recognition site starts at (position - cutTop)
                // For EcoRI at cut position 101: GAATTC is at 0-indexed 100-105
                // cutTop=1, so site starts at 101 - 1 = 100
                let topCutPos = position  // Position where top strand cuts (0-indexed)
                
                // Show context bases ending right before the top cut
                let contextStart = max(0, topCutPos - contextBases)
                let contextEnd = topCutPos  // exclusive end
                let leftContext = extractSequence(start: contextStart, end: contextEnd)
                let leftContextComp = complementSequence(leftContext)
                
                // Overhang: bases between top cut and bottom cut (on the bottom strand)
                let overhangStart = position
                let overhangEnd = siteStart + cutBottom
                let overhang = extractSequence(start: overhangStart, end: overhangEnd)
                let overhangComp = complementSequence(overhang)
                
                // Top strand: just context, ends at cut
                // Bottom strand: context + overhang complement (continues past the top cut)
                topStrand = leftContext + String(repeating: " ", count: overhangSize) + " -5'"
                bottomStrand = leftContextComp + overhangComp + " -3'"
            } else {
                // Right side: top strand starts before bottom strand
                // Top strand begins at its cut point, bottom strand is recessed
                
                // The overhang is what the top strand has before the bottom strand starts
                let overhangStart = position
                let overhangEnd = siteStart + cutBottom
                let overhang = extractSequence(start: overhangStart, end: overhangEnd)
                
                let contextStart = siteStart + cutBottom
                let contextEnd = min(sequence.count, siteStart + cutBottom + contextBases)
                let rightContext = extractSequence(start: contextStart, end: contextEnd)
                let rightContextComp = complementSequence(rightContext)
                
                // Top strand: overhang + context
                // Bottom strand: spaces to align + context complement
                topStrand = overhang + rightContext + " -5'"
                bottomStrand = String(repeating: " ", count: overhangSize) + rightContextComp + " -3'"
            }
        } else {
            // 3' overhang - cutTop > cutBottom
            // bottom cuts first, top cuts later
            // For 3' overhang: bottom strand is shorter on left, top is shorter on right
            
            if isLeftEnd {
                // Left side: bottom strand is shorter, top strand extends with overhang
                let bottomCutPos = siteStart + cutBottom  // Bottom strand cuts first
                
                // Show context bases ending at the bottom cut
                let contextStart = max(0, bottomCutPos - contextBases)
                let contextEnd = bottomCutPos
                let leftContext = extractSequence(start: contextStart, end: contextEnd)
                let leftContextComp = complementSequence(leftContext)
                
                // Overhang: bases between bottom cut and top cut (on the top strand)
                let overhangStart = siteStart + cutBottom
                let overhangEnd = siteStart + cutTop
                let overhang = extractSequence(start: overhangStart, end: overhangEnd)
                
                // Top strand: context + overhang (extends further)
                // Bottom strand: context complement + spaces
                topStrand = leftContext + overhang + " -5'"
                bottomStrand = leftContextComp + String(repeating: " ", count: overhangSize) + " -3'"
            } else {
                // Right side: top strand is shorter, bottom strand has overhang
                
                // The overhang is what the bottom strand has before the top strand starts
                let overhangStart = siteStart + cutBottom
                let overhangEnd = position  // Top cuts at position
                let overhang = extractSequence(start: overhangStart, end: overhangEnd)
                let overhangComp = complementSequence(overhang)
                
                let contextStart = position
                let contextEnd = min(sequence.count, position + contextBases)
                let rightContext = extractSequence(start: contextStart, end: contextEnd)
                let rightContextComp = complementSequence(rightContext)
                
                // Top strand: spaces + context (recessed)
                // Bottom strand: overhang complement + context complement (extends left)
                topStrand = String(repeating: " ", count: overhangSize) + rightContext + " -5'"
                bottomStrand = overhangComp + rightContextComp + " -3'"
            }
        }
        
        return DNAVisualization(topStrand: topStrand, bottomStrand: bottomStrand)
    }
    
    private func extractSequence(start: Int, end: Int) -> String {
        guard start >= 0 && end <= sequence.count && start < end else {
            return ""
        }
        let startIdx = sequence.index(sequence.startIndex, offsetBy: start)
        let endIdx = sequence.index(sequence.startIndex, offsetBy: end)
        return String(sequence[startIdx..<endIdx])
    }
}

// MARK: - Preview
#Preview {
    NavigationStack {
        CutSiteDetailView(
            cutSite: CutSite(
                enzyme: Enzyme(
                    name: "EcoRI",
                    site: "GAATTC",
                    cutIndexTop: 1,
                    cutIndexBottom: 5,
                    overhangType: .fivePrime,
                    notes: nil
                ),
                position: 101,
                strand: .top
            ),
            index: 0,
            sequence: "ATGCATGCATGCGAATTCATGCATGCATGCATGCATGCATGC"
        )
    }
}


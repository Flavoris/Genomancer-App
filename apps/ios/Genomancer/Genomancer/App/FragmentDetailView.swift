import SwiftUI
import DigestCore

struct FragmentDetailView: View {
    let fragment: Fragment
    let index: Int
    
    @State private var showFullSequence = false
    
    var body: some View {
        ScrollView {
            VStack(alignment: .leading, spacing: 20) {
                // Fragment Summary Card
                fragmentSummaryCard
                
                // Ends Information
                endsInformationCard
                
                // Sequence Display
                sequenceCard
            }
            .padding()
        }
        .background(Color.genomancerBackground)
        .navigationTitle("Fragment \(index + 1)")
        .navigationBarTitleDisplayMode(.inline)
        .toolbarColorScheme(.dark, for: .navigationBar)
    }
    
    private var fragmentSummaryCard: some View {
        VStack(alignment: .leading, spacing: 12) {
            Text("Fragment Information")
                .font(.headline)
                .foregroundColor(.primary)
            
            Divider()
            
            HStack {
                Text("Length:")
                    .foregroundColor(.secondary)
                Spacer()
                Text("\(fragment.length) bp")
                    .font(.headline)
                    .foregroundColor(.genomancerAccent)
            }
            
            HStack {
                Text("Position:")
                    .foregroundColor(.secondary)
                Spacer()
                Text("[\(fragment.start) .. \(fragment.end))")
                    .font(.system(.body, design: .monospaced))
                    .foregroundColor(.primary)
            }
        }
        .padding()
        .background(Color(.systemBackground))
        .cornerRadius(12)
    }
    
    private var endsInformationCard: some View {
        VStack(alignment: .leading, spacing: 16) {
            Text("Fragment Ends")
                .font(.headline)
                .foregroundColor(.primary)
            
            Divider()
            
            // Left End (5')
            VStack(alignment: .leading, spacing: 8) {
                HStack {
                    Text("5' End (Left)")
                        .font(.subheadline)
                        .fontWeight(.semibold)
                        .foregroundColor(.genomancerAccent)
                    Spacer()
                }
                
                endDetailsWithVisualization(fragment.leftEnd, isLeftEnd: true)
            }
            .padding()
            .background(Color(.systemGray6))
            .cornerRadius(8)
            
            // Right End (3')
            VStack(alignment: .leading, spacing: 8) {
                HStack {
                    Text("3' End (Right)")
                        .font(.subheadline)
                        .fontWeight(.semibold)
                        .foregroundColor(.genomancerAccent)
                    Spacer()
                }
                
                endDetailsWithVisualization(fragment.rightEnd, isLeftEnd: false)
            }
            .padding()
            .background(Color(.systemGray6))
            .cornerRadius(8)
        }
        .padding()
        .background(Color(.systemBackground))
        .cornerRadius(12)
    }
    
    private func endDetailsWithVisualization(_ end: EndInfo, isLeftEnd: Bool) -> some View {
        VStack(alignment: .leading, spacing: 8) {
            if let enzyme = end.sourceEnzyme {
                HStack {
                    Text("Enzyme:")
                        .foregroundColor(.secondary)
                        .font(.caption)
                    Text(enzyme)
                        .font(.caption)
                        .fontWeight(.medium)
                        .foregroundColor(.primary)
                }
            }
            
            HStack {
                Text("Type:")
                    .foregroundColor(.secondary)
                    .font(.caption)
                Text(overhangDescription(end.overhangType))
                    .font(.caption)
                    .foregroundColor(.primary)
            }
            
            if let seq = end.overhangSeq5to3, !seq.isEmpty {
                VStack(alignment: .leading, spacing: 4) {
                    Text("Overhang:")
                        .foregroundColor(.secondary)
                        .font(.caption)
                    
                    // DNA visualization
                    if let visualization = createDNAVisualization(end: end, isLeftEnd: isLeftEnd) {
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
    }
    
    private var sequenceCard: some View {
        VStack(alignment: .leading, spacing: 12) {
            HStack {
                Text("DNA Sequence")
                    .font(.headline)
                    .foregroundColor(.primary)
                
                Spacer()
                
                if let seq = fragment.sequence, seq.count > 120 {
                    Button(action: {
                        showFullSequence.toggle()
                    }) {
                        Text(showFullSequence ? "Show Less" : "Show All")
                            .font(.caption)
                            .foregroundColor(.genomancerAccent)
                    }
                }
            }
            
            Divider()
            
            if let seq = fragment.sequence {
                sequenceDisplay(seq)
            } else {
                Text("Sequence not available")
                    .foregroundColor(.secondary)
                    .font(.caption)
                    .italic()
            }
            
            // Copy button
            if let seq = fragment.sequence {
                Button(action: {
                    UIPasteboard.general.string = seq
                }) {
                    HStack {
                        Image(systemName: "doc.on.doc")
                        Text("Copy Sequence")
                    }
                    .font(.caption)
                    .foregroundColor(.genomancerAccent)
                }
                .buttonStyle(.bordered)
            }
        }
        .padding()
        .background(Color(.systemBackground))
        .cornerRadius(12)
    }
    
    private func sequenceDisplay(_ sequence: String) -> some View {
        VStack(alignment: .leading, spacing: 4) {
            if showFullSequence || sequence.count <= 120 {
                // Show full formatted sequence
                Text(formatSequence(sequence, lineLength: 60, groupSize: 10))
                    .font(.system(.caption, design: .monospaced))
                    .foregroundColor(.primary)
                    .textSelection(.enabled)
            } else {
                // Show truncated sequence
                let preview = sequence.prefix(60).description + " ... " + sequence.suffix(60).description
                Text(formatSequence(preview, lineLength: 60, groupSize: 10))
                    .font(.system(.caption, design: .monospaced))
                    .foregroundColor(.primary)
                    .textSelection(.enabled)
                
                Text("(\(sequence.count - 120) more bases hidden)")
                    .font(.caption2)
                    .foregroundColor(.secondary)
                    .italic()
            }
        }
    }
    
    private func formatSequence(_ seq: String, lineLength: Int, groupSize: Int) -> String {
        var result = ""
        var position = 0
        
        while position < seq.count {
            let lineEnd = min(position + lineLength, seq.count)
            let lineStart = seq.index(seq.startIndex, offsetBy: position)
            let lineEndIndex = seq.index(seq.startIndex, offsetBy: lineEnd)
            let line = String(seq[lineStart..<lineEndIndex])
            
            // Add grouped bases
            var groupPos = 0
            for char in line {
                if groupPos > 0 && groupPos % groupSize == 0 {
                    result += " "
                }
                result.append(char)
                groupPos += 1
            }
            
            result += "\n"
            position = lineEnd
        }
        
        return result
    }
    
    private func overhangDescription(_ type: OverhangType) -> String {
        switch type {
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
    
    // MARK: - DNA Visualization Helpers
    
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
    
    private func createDNAVisualization(end: EndInfo, isLeftEnd: Bool) -> DNAVisualization? {
        guard let overhangSeq = end.overhangSeq5to3,
              !overhangSeq.isEmpty,
              let sequence = fragment.sequence else {
            return nil
        }
        
        let contextBases = 3
        let overhangLen = overhangSeq.count
        
        var contextSeq: String
        
        if isLeftEnd {
            // For left end, get first few bases from fragment
            let endIndex = min(contextBases, sequence.count)
            contextSeq = String(sequence.prefix(endIndex))
        } else {
            // For right end, get last few bases from fragment
            contextSeq = String(sequence.suffix(contextBases))
        }
        
        let contextComplement = complementSequence(contextSeq)
        let overhangComplement = complementSequence(overhangSeq)
        
        var topStrand: String
        var bottomStrand: String
        
        switch end.overhangType {
        case .fivePrime:
            // 5' overhang extends on top strand
            if isLeftEnd {
                topStrand = "5'- \(overhangSeq) \(contextSeq)"
                bottomStrand = "3'-" + String(repeating: " ", count: overhangLen + 1) + contextComplement
            } else {
                topStrand = "   \(contextSeq) \(overhangSeq) -5'"
                bottomStrand = "   " + contextComplement + String(repeating: " ", count: overhangLen + 1) + "-3'"
            }
            
        case .threePrime:
            // 3' overhang extends on bottom strand
            if isLeftEnd {
                topStrand = "5'- " + contextSeq
                bottomStrand = "3'- \(overhangComplement) \(contextComplement)"
            } else {
                topStrand = "   " + contextSeq + String(repeating: " ", count: overhangLen + 1) + "-5'"
                bottomStrand = "   \(contextComplement) \(overhangComplement) -3'"
            }
            
        case .blunt:
            // Blunt end - both strands same length
            if isLeftEnd {
                topStrand = "5'- \(contextSeq)"
                bottomStrand = "3'- \(contextComplement)"
            } else {
                topStrand = "   \(contextSeq) -5'"
                bottomStrand = "   \(contextComplement) -3'"
            }
            
        case .unknown:
            return nil
        }
        
        return DNAVisualization(topStrand: topStrand, bottomStrand: bottomStrand)
    }
}

// MARK: - Preview
#Preview {
    NavigationStack {
        FragmentDetailView(
            fragment: Fragment(
                start: 0,
                end: 101,
                length: 101,
                leftEnd: EndInfo(
                    overhangType: .blunt,
                    overhangSeq5to3: nil,
                    sourceEnzyme: nil
                ),
                rightEnd: EndInfo(
                    overhangType: .fivePrime,
                    overhangSeq5to3: "AATT",
                    sourceEnzyme: "EcoRI"
                ),
                sequence: "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
            ),
            index: 0
        )
    }
}


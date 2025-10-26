import SwiftUI
import DigestCore

struct FragmentList: View {
    let fragments: [Fragment]
    let fullSequence: String
    
    var body: some View {
        List(Array(fragments.enumerated()), id: \.element) { index, f in
            NavigationLink(destination: FragmentDetailView(fragment: f, index: index)) {
                FragmentRowView(fragment: f, index: index)
            }
        }
        .scrollContentBackground(.hidden)
        .background(Color.genomancerBackground)
        .navigationTitle("Fragments")
        .toolbarColorScheme(.dark, for: .navigationBar)
    }
}

struct FragmentRowView: View {
    let fragment: Fragment
    let index: Int
    
    var body: some View {
        VStack(alignment: .leading, spacing: 8) {
            // Fragment label and position
            HStack {
                Text("Fragment \(index + 1)")
                    .font(.subheadline)
                    .fontWeight(.semibold)
                    .foregroundColor(.genomancerAccent)
                
                Spacer()
                
                Text("[\(fragment.start) .. \(fragment.end))")
                    .font(.caption)
                    .foregroundColor(.genomancerSecondaryText)
                    .dynamicTypeSize(...DynamicTypeSize.accessibility3)
            }
            
            // Fragment length
            HStack {
                Text("\(fragment.length) bp")
                    .font(.headline)
                    .foregroundColor(.genomancerText)
                    .dynamicTypeSize(...DynamicTypeSize.accessibility3)
                
                Spacer()
            }
            
            // Overhang info
            HStack(spacing: 12) {
                // Left end
                if let enzyme = fragment.leftEnd.sourceEnzyme {
                    HStack(spacing: 4) {
                        Text("5':")
                            .font(.caption2)
                            .foregroundColor(.genomancerSecondaryText)
                        Text(enzyme)
                            .font(.caption2)
                            .foregroundColor(.genomancerAccent)
                        if let seq = fragment.leftEnd.overhangSeq5to3 {
                            Text("(\(seq))")
                                .font(.system(.caption2, design: .monospaced))
                                .foregroundColor(.genomancerSecondaryText)
                        }
                    }
                }
                
                // Right end
                if let enzyme = fragment.rightEnd.sourceEnzyme {
                    HStack(spacing: 4) {
                        Text("3':")
                            .font(.caption2)
                            .foregroundColor(.genomancerSecondaryText)
                        Text(enzyme)
                            .font(.caption2)
                            .foregroundColor(.genomancerAccent)
                        if let seq = fragment.rightEnd.overhangSeq5to3 {
                            Text("(\(seq))")
                                .font(.system(.caption2, design: .monospaced))
                                .foregroundColor(.genomancerSecondaryText)
                        }
                    }
                }
            }
            
            // Sequence preview
            if let seq = fragment.sequence {
                SequencePreview(sequence: seq)
            }
        }
        .padding(.vertical, 4)
        .accessibilityElement(children: .combine)
        .accessibilityLabel("Fragment \(index + 1)")
        .accessibilityValue("\(fragment.length) base pairs, from position \(fragment.start) to \(fragment.end)")
        .accessibilityHint("Left overhang: \(overhangDescription(fragment.leftEnd)), Right overhang: \(overhangDescription(fragment.rightEnd)). Tap for details.")
    }
    
    private func overhangDescription(_ end: EndInfo) -> String {
        switch end.overhangType {
        case .blunt:
            return "blunt"
        case .fivePrime:
            return "5' overhang" + (end.overhangSeq5to3 != nil ? " \(end.overhangSeq5to3!)" : "")
        case .threePrime:
            return "3' overhang" + (end.overhangSeq5to3 != nil ? " \(end.overhangSeq5to3!)" : "")
        case .unknown:
            return "unknown"
        }
    }
}

struct SequencePreview: View {
    let sequence: String
    
    private var previewText: String {
        if sequence.count <= 40 {
            return sequence
        } else {
            return sequence.prefix(20).description + "..." + sequence.suffix(17).description
        }
    }
    
    var body: some View {
        HStack(spacing: 4) {
            Image(systemName: "doc.text")
                .font(.caption2)
                .foregroundColor(.genomancerSecondaryText)
            
            Text(previewText)
                .font(.system(.caption2, design: .monospaced))
                .foregroundColor(.genomancerSecondaryText)
                .lineLimit(1)
        }
        .padding(.horizontal, 8)
        .padding(.vertical, 4)
        .background(Color.genomancerSecondaryBackground.opacity(0.3))
        .cornerRadius(6)
    }
}


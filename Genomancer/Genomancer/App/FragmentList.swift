import SwiftUI
import DigestCore

struct FragmentList: View {
    let fragments: [Fragment]
    var body: some View {
        List(Array(fragments.enumerated()), id: \.element) { index, f in
            HStack {
                Text("\(f.length) bp")
                    .font(.headline)
                    .dynamicTypeSize(...DynamicTypeSize.accessibility3)
                Spacer()
                Text("[\(f.start) .. \(f.end))")
                    .font(.caption)
                        .foregroundColor(.genomancerSecondaryText)
                    .dynamicTypeSize(...DynamicTypeSize.accessibility3)
            }
            .accessibilityElement(children: .combine)
            .accessibilityLabel("Fragment \(index + 1)")
            .accessibilityValue("\(f.length) base pairs, from position \(f.start) to \(f.end)")
            .accessibilityHint("Left overhang: \(overhangDescription(f.leftEnd)), Right overhang: \(overhangDescription(f.rightEnd))")
        }
        .scrollContentBackground(.hidden)
        .background(Color.genomancerBackground)
        .navigationTitle("Fragments")
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


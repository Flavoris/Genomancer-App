import SwiftUI
import DigestCore

struct GelView: View {
    let fragments: [Fragment]
    
    @State private var gelText: String = ""
    
    var body: some View {
        ScrollView(.vertical) {
            Text(gelText)
                .font(.system(.caption, design: .monospaced))
                .foregroundColor(.genomancerText)
                .dynamicTypeSize(...DynamicTypeSize.accessibility1)
                .minimumScaleFactor(0.5)
                .lineLimit(nil)
                .padding()
                .accessibilityLabel("Gel electrophoresis visualization")
                .accessibilityValue(gelAccessibilityDescription)
        }
        .background(Color.genomancerBackground)
        .navigationTitle("Gel Electrophoresis")
        .navigationBarTitleDisplayMode(.inline)
        .onAppear(perform: generateGel)
    }
    
    private var gelAccessibilityDescription: String {
        let lengths = fragments.map { $0.length }.sorted(by: >)
        if lengths.isEmpty {
            return "No fragments"
        }
        let fragmentDesc = lengths.prefix(5).map { "\($0) base pairs" }.joined(separator: ", ")
        let suffix = lengths.count > 5 ? " and \(lengths.count - 5) more" : ""
        return "Showing \(lengths.count) fragments: \(fragmentDesc)\(suffix)"
    }
    
    private func generateGel() {
        let lengths = fragments.map { $0.length }
        gelText = asciiGel(fragments: lengths)
    }
}


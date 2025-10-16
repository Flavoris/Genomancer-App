import SwiftUI
import DigestCore

struct GelView: View {
    let fragments: [Fragment]
    
    @State private var gelText: String = ""
    
    var body: some View {
        ScrollView(.vertical) {
            Text(gelText)
                .font(.system(.caption, design: .monospaced))
                .minimumScaleFactor(0.5)
                .lineLimit(nil)
                .padding()
        }
        .navigationTitle("Gel Electrophoresis")
        .navigationBarTitleDisplayMode(.inline)
        .onAppear(perform: generateGel)
    }
    
    private func generateGel() {
        let lengths = fragments.map { $0.length }
        gelText = asciiGel(fragments: lengths)
    }
}


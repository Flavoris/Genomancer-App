import SwiftUI
import DigestCore

struct GelView: View {
    let fragments: [Fragment]
    
    var body: some View {
        GelRenderView(fragments: fragments, style: .photo)
            .background(Color.genomancerBackground)
            .navigationTitle("Gel Electrophoresis")
            .navigationBarTitleDisplayMode(.inline)
            .toolbarColorScheme(.dark, for: .navigationBar)
    }
}


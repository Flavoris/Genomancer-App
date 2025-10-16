import SwiftUI
import DigestCore

struct FragmentList: View {
    let fragments: [Fragment]
    var body: some View {
        List(fragments, id: \.self) { f in
            HStack {
                Text("\(f.length) bp").font(.headline)
                Spacer()
                Text("[\(f.start) .. \(f.end))").font(.caption).foregroundStyle(.secondary)
            }
        }
        .navigationTitle("Fragments")
    }
}


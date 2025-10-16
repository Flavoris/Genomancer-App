import SwiftUI
import DigestCore

struct HomeView: View {
    @State private var sequence: String = ""
    @State private var allEnzymes: [Enzyme] = []
    @State private var selected: Set<Enzyme> = []
    @State private var circular = false
    @State private var results: [Fragment] = []

    var body: some View {
        NavigationStack {
            Form {
                Section("Sequence (FASTA or raw)") {
                    TextEditor(text: $sequence)
                        .frame(minHeight: 140)
                        .font(.system(.body, design: .monospaced))
                }
                Section("Options") {
                    Toggle("Circular (plasmid)", isOn: $circular)
                    NavigationLink("Choose Enzymes") { EnzymePicker(all: allEnzymes, selected: $selected) }
                }
                Button("Digest") { runDigest() }.buttonStyle(.borderedProminent)
                if !results.isEmpty {
                    NavigationLink("View Fragments") { FragmentList(fragments: results) }
                }
            }
            .navigationTitle("Genomancer")
            .onAppear(perform: loadEnzymes)
        }
    }

    func loadEnzymes() {
        guard let url = Bundle.main.url(forResource: "enzymes", withExtension: "json") else { return }
        if let loaded = try? EnzymeLoader.load(from: url) { self.allEnzymes = loaded }
    }

    func runDigest() {
        let dna = parseFASTA(sequence)
        let engine = DigestEngine(sequence: dna, enzymes: Array(selected))
        results = engine.digest(options: .init(circular: circular, returnSequences: false))
    }

    func parseFASTA(_ s: String) -> String {
        s.split(separator: "\n").filter{ !$0.hasPrefix(">") }
            .joined().replacingOccurrences(of: "\\s", with: "", options: .regularExpression)
    }
}


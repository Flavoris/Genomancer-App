import SwiftUI
import DigestCore
import UniformTypeIdentifiers

struct HomeView: View {
    @State private var sequence: String = ""
    @State private var allEnzymes: [Enzyme] = []
    @State private var selected: Set<Enzyme> = []
    @State private var circular = false
    @State private var results: [Fragment] = []
    @State private var csvExportData: ExportData?
    @State private var genbankExportData: ExportData?

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
                    NavigationLink("View Map") { 
                        MapView(sequence: parseFASTA(sequence), enzymes: Array(selected), circular: circular) 
                    }
                    NavigationLink("View Gel") { GelView(fragments: results) }
                    
                    Section("Export") {
                        if let csvData = csvExportData {
                            ShareLink(
                                item: csvData.data,
                                preview: SharePreview(csvData.filename, image: Image(systemName: "doc.text"))
                            ) {
                                Label("Export CSV", systemImage: "tablecells")
                            }
                        }
                        
                        if let gbData = genbankExportData {
                            ShareLink(
                                item: gbData.data,
                                preview: SharePreview(gbData.filename, image: Image(systemName: "doc.text"))
                            ) {
                                Label("Export GenBank", systemImage: "doc.text")
                            }
                        }
                    }
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
        
        // Generate export data
        prepareExports(dna: dna)
    }

    func parseFASTA(_ s: String) -> String {
        s.split(separator: "\n").filter{ !$0.hasPrefix(">") }
            .joined().replacingOccurrences(of: "\\s", with: "", options: .regularExpression)
    }
    
    func prepareExports(dna: String) {
        // Extract locus name from FASTA header or use timestamp
        let locusName = extractLocusName()
        let timestamp = formattedTimestamp()
        let filename = "\(locusName)_\(timestamp)"
        
        // Generate CSV
        let csvContent = exportCSV(fragments: results)
        csvExportData = ExportData(
            filename: "\(filename)_fragments.csv",
            data: csvContent,
            utType: .commaSeparatedText
        )
        
        // Generate GenBank
        let features = results.enumerated().map { (index, frag) -> (Range<Int>, String) in
            let label = "fragment_\(index + 1)"
            return (frag.start..<frag.end, label)
        }
        let gbContent = exportGenBank(sequence: dna, locus: locusName, features: features)
        genbankExportData = ExportData(
            filename: "\(filename).gb",
            data: gbContent,
            utType: .plainText
        )
    }
    
    func extractLocusName() -> String {
        // Try to extract name from FASTA header
        let lines = sequence.split(separator: "\n")
        if let firstLine = lines.first, firstLine.hasPrefix(">") {
            let header = String(firstLine.dropFirst())
            let cleaned = header.split(separator: " ").first ?? "sequence"
            return String(cleaned).replacingOccurrences(of: "[^A-Za-z0-9_]", with: "_", options: .regularExpression)
        }
        return "sequence"
    }
    
    func formattedTimestamp() -> String {
        let formatter = DateFormatter()
        formatter.dateFormat = "yyyyMMdd_HHmmss"
        return formatter.string(from: Date())
    }
}

struct ExportData: Identifiable, Transferable {
    let id = UUID()
    let filename: String
    let data: String
    let utType: UTType
    
    static var transferRepresentation: some TransferRepresentation {
        DataRepresentation(exportedContentType: .plainText) { exportData in
            Data(exportData.data.utf8)
        }
    }
    
    var suggestedName: String { filename }
}


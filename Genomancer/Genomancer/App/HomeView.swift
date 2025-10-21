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
    @State private var showValidationAlert = false
    @State private var validationMessage = ""
    @State private var showCopyAlert = false
    @State private var copyMessage = ""
    @State private var showFileImporter = false
    @State private var showFileImportError = false
    @State private var fileImportErrorMessage = ""

    var body: some View {
        NavigationStack {
            ScrollView {
                VStack(spacing: 0) {
                    logoHeader
                    formContent
                }
            }
            .background(Color.genomancerBackground)
            .navigationTitle("")
            .navigationBarTitleDisplayMode(.inline)
            .onAppear(perform: loadEnzymes)
            .alert("Validation Warning", isPresented: $showValidationAlert) {
                Button("OK", role: .cancel) { }
            } message: {
                Text(validationMessage)
            }
            .alert("Copied", isPresented: $showCopyAlert) {
                Button("OK", role: .cancel) { }
            } message: {
                Text(copyMessage)
            }
            .alert("Import Error", isPresented: $showFileImportError) {
                Button("OK", role: .cancel) { }
            } message: {
                Text(fileImportErrorMessage)
            }
            .fileImporter(
                isPresented: $showFileImporter,
                allowedContentTypes: [.plainText, .text, UTType(filenameExtension: "fasta") ?? .plainText, UTType(filenameExtension: "fa") ?? .plainText, UTType(filenameExtension: "fna") ?? .plainText],
                allowsMultipleSelection: false
            ) { result in
                handleFileImport(result: result)
            }
        }
    }
    
    private var logoHeader: some View {
        VStack(spacing: -8) {
            Image("Logo")
                .resizable()
                .aspectRatio(contentMode: .fit)
                .frame(width: 100, height: 100)
                .accessibilityLabel("Genomancer logo - wizard hat with DNA helix")
            
            Image("genomancer-font")
                .resizable()
                .aspectRatio(contentMode: .fit)
                .frame(height: 160)
                .padding(.horizontal, 10)
                .accessibilityLabel("Genomancer DNA Restriction Analysis")
        }
        .frame(maxWidth: .infinity)
        .padding(.horizontal)
        .background(Color.genomancerBackground)
    }
    
    private var formContent: some View {
        VStack(spacing: 20) {
            // DNA Sequence Section
            VStack(alignment: .leading, spacing: 12) {
                Text("DNA Sequence:")
                    .foregroundColor(.genomancerText)
                    .font(.headline)
                    .padding(.horizontal)
                
                VStack(spacing: 12) {
                    TextEditor(text: $sequence)
                        .frame(minHeight: 140)
                        .font(.system(.body, design: .monospaced))
                        .dynamicTypeSize(...DynamicTypeSize.accessibility2)
                        .onChange(of: sequence) { newValue in
                            validateSequence(newValue)
                        }
                        .accessibilityLabel("DNA sequence input")
                        .accessibilityHint("Enter DNA sequence in FASTA or raw format")
                        .padding(8)
                        .background(Color(.systemBackground))
                        .cornerRadius(8)
                    
                    if isFASTAFormat {
                        Label("FASTA format detected", systemImage: "doc.text")
                            .font(.caption)
                            .foregroundColor(.genomancerSecondaryText)
                            .accessibilityLabel("FASTA format detected in sequence")
                            .frame(maxWidth: .infinity, alignment: .leading)
                    }
                    
                    Button(action: { showFileImporter = true }) {
                        Label("Import FASTA File", systemImage: "doc.badge.plus")
                            .frame(maxWidth: .infinity)
                    }
                    .buttonStyle(.bordered)
                    .accessibilityLabel("Import FASTA file")
                    .accessibilityHint("Opens file picker to select a FASTA file")
                }
                .padding(.horizontal)
            }
            
            // Enzyme Selection Section
            VStack(spacing: 12) {
                HStack(spacing: 8) {
                    Text("Circular (plasmid)")
                    Spacer()
                    Toggle("", isOn: $circular)
                        .labelsHidden()
                }
                
                
                NavigationLink {
                    EnzymePicker(all: allEnzymes, selected: $selected)
                } label: {
                    HStack {
                        Text("Choose Enzyme(s)")
                            .foregroundColor(.primary)
                        Spacer()
                        Image(systemName: "chevron.right")
                            .foregroundColor(.secondary)
                            .font(.system(size: 14, weight: .semibold))
                    }
                }
            }
            .padding()
            .background(Color(.systemBackground))
            .cornerRadius(8)
            .padding(.horizontal)
            
            // Digest Button
            Button("Digest") { runDigest() }
                .buttonStyle(GenomancerProminentButtonStyle())
                .disabled(!canDigest)
                .frame(maxWidth: .infinity)
                .padding(.horizontal)
            
            // Results Section
            if !results.isEmpty {
                VStack(spacing: 12) {
                    NavigationLink("View Fragments") { 
                        FragmentList(fragments: results) 
                    }
                    .padding()
                    .background(Color(.systemBackground))
                    .cornerRadius(8)
                    
                    NavigationLink("View Map") { 
                        MapView(sequence: parseFASTA(sequence), enzymes: Array(selected), circular: circular) 
                    }
                    .padding()
                    .background(Color(.systemBackground))
                    .cornerRadius(8)
                    
                    NavigationLink("View Gel") { 
                        GelView(fragments: results) 
                    }
                    .padding()
                    .background(Color(.systemBackground))
                    .cornerRadius(8)
                }
                .padding(.horizontal)
                
                // Export Section
                VStack(alignment: .leading, spacing: 12) {
                    Text("Export")
                        .foregroundColor(.genomancerText)
                        .font(.headline)
                        .padding(.horizontal)
                    
                    VStack(spacing: 12) {
                        if let csvData = csvExportData {
                            HStack {
                                ShareLink(
                                    item: csvData.data,
                                    preview: SharePreview(csvData.filename, image: Image(systemName: "doc.text"))
                                ) {
                                    Label("Export CSV", systemImage: "tablecells")
                                }
                                Spacer()
                                Button(action: { copyToClipboard(csvData.data, format: "CSV") }) {
                                    Image(systemName: "doc.on.doc")
                                        .foregroundColor(.genomancerAccent)
                                }
                                .buttonStyle(.borderless)
                            }
                            .padding()
                            .background(Color(.systemBackground))
                            .cornerRadius(8)
                        }
                        
                        if let gbData = genbankExportData {
                            HStack {
                                ShareLink(
                                    item: gbData.data,
                                    preview: SharePreview(gbData.filename, image: Image(systemName: "doc.text"))
                                ) {
                                    Label("Export GenBank", systemImage: "doc.text")
                                }
                                Spacer()
                                Button(action: { copyToClipboard(gbData.data, format: "GenBank") }) {
                                    Image(systemName: "doc.on.doc")
                                        .foregroundColor(.genomancerAccent)
                                }
                                .buttonStyle(.borderless)
                            }
                            .padding()
                            .background(Color(.systemBackground))
                            .cornerRadius(8)
                        }
                    }
                    .padding(.horizontal)
                }
            }
        }
        .padding(.vertical)
    }
    
    var isFASTAFormat: Bool {
        sequence.trimmingCharacters(in: .whitespacesAndNewlines).hasPrefix(">")
    }
    
    var canDigest: Bool {
        let dna = parseFASTA(sequence)
        return !dna.isEmpty && !selected.isEmpty
    }

    func loadEnzymes() {
        guard let url = Bundle.main.url(forResource: "enzymes", withExtension: "json") else { return }
        if let loaded = try? EnzymeLoader.load(from: url) { self.allEnzymes = loaded }
    }

    func runDigest() {
        let dna = parseFASTA(sequence)
        
        // Guard against empty input
        guard !dna.isEmpty else {
            validationMessage = "Please enter a DNA sequence."
            showValidationAlert = true
            return
        }
        
        guard !selected.isEmpty else {
            validationMessage = "Please select at least one restriction enzyme."
            showValidationAlert = true
            return
        }
        
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
    
    func validateSequence(_ input: String) {
        let dna = parseFASTA(input)
        guard !dna.isEmpty else { return }
        
        let validChars = CharacterSet(charactersIn: "ATCGNRYSWKMBDHVatcgryswkmbdhv")
        let inputChars = CharacterSet(charactersIn: dna)
        
        if !inputChars.isSubset(of: validChars) {
            validationMessage = "Warning: Sequence contains illegal characters. Only IUPAC DNA codes (A,T,C,G,N,R,Y,S,W,K,M,B,D,H,V) are valid."
            showValidationAlert = true
        }
    }
    
    func copyToClipboard(_ content: String, format: String) {
        UIPasteboard.general.string = content
        copyMessage = "\(format) data copied to clipboard"
        showCopyAlert = true
    }
    
    func handleFileImport(result: Result<[URL], Error>) {
        switch result {
        case .success(let urls):
            guard let selectedFile = urls.first else {
                fileImportErrorMessage = "No file was selected."
                showFileImportError = true
                return
            }
            
            // Start accessing the security-scoped resource
            guard selectedFile.startAccessingSecurityScopedResource() else {
                fileImportErrorMessage = "Unable to access the selected file."
                showFileImportError = true
                return
            }
            
            defer { selectedFile.stopAccessingSecurityScopedResource() }
            
            do {
                let fileContent = try String(contentsOf: selectedFile, encoding: .utf8)
                
                // Validate that it's a valid FASTA or plain sequence
                let trimmedContent = fileContent.trimmingCharacters(in: .whitespacesAndNewlines)
                if trimmedContent.isEmpty {
                    fileImportErrorMessage = "The selected file is empty."
                    showFileImportError = true
                    return
                }
                
                // Update the sequence field with the file content
                sequence = fileContent
                
            } catch {
                fileImportErrorMessage = "Failed to read file: \(error.localizedDescription)"
                showFileImportError = true
            }
            
        case .failure(let error):
            fileImportErrorMessage = "Failed to import file: \(error.localizedDescription)"
            showFileImportError = true
        }
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


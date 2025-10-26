import SwiftUI
import DigestCore
import UniformTypeIdentifiers

struct HomeView: View {
    @State private var sequence: String = ""
    @State private var allEnzymes: [Enzyme] = []
    @State private var selected: Set<Enzyme> = []
    @State private var circular = false
    @State private var results: [Fragment] = []
    @State private var showValidationAlert = false
    @State private var validationMessage = ""
    @State private var showFileImporter = false
    @State private var showFileImportError = false
    @State private var fileImportErrorMessage = ""
    @FocusState private var isSequenceFieldFocused: Bool

    var body: some View {
        NavigationStack {
            ScrollView {
                VStack(spacing: 0) {
                    logoHeader
                    formContent
                }
            }
            .background(Color.genomancerBackground)
            .onTapGesture {
                isSequenceFieldFocused = false
            }
            .navigationTitle("")
            .navigationBarTitleDisplayMode(.inline)
            .onAppear(perform: loadEnzymes)
            .alert("Validation Warning", isPresented: $showValidationAlert) {
                Button("OK", role: .cancel) { }
            } message: {
                Text(validationMessage)
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
                        .frame(height: 120)
                        .font(.system(.body, design: .monospaced))
                        .dynamicTypeSize(...DynamicTypeSize.accessibility2)
                        .focused($isSequenceFieldFocused)
                        .onChange(of: sequence) { newValue in
                            validateSequence(newValue)
                        }
                        .accessibilityLabel("DNA sequence input")
                        .accessibilityHint("Enter DNA sequence in FASTA or raw format. Scrollable for longer sequences.")
                        .padding(8)
                        .background(Color(.systemBackground))
                        .cornerRadius(8)
                        .toolbar {
                            ToolbarItemGroup(placement: .keyboard) {
                                Spacer()
                                Button("Done") {
                                    isSequenceFieldFocused = false
                                }
                            }
                        }
                    
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
                
                Divider()
                
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
                
                // Display selected enzymes
                if !selected.isEmpty {
                    selectedEnzymesView
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
                LazyVGrid(columns: [
                    GridItem(.flexible(), spacing: 12),
                    GridItem(.flexible(), spacing: 12)
                ], spacing: 12) {
                    NavigationLink("View Fragments") { 
                        FragmentList(fragments: results, fullSequence: parseFASTA(sequence)) 
                    }
                    .padding()
                    .background(Color(.systemBackground))
                    .cornerRadius(8)
                    
                    if circular {
                        NavigationLink("View Map") { 
                            MapView(sequence: parseFASTA(sequence), enzymes: Array(selected), circular: circular) 
                        }
                        .padding()
                        .background(Color(.systemBackground))
                        .cornerRadius(8)
                    }
                    
                    NavigationLink("View Gel") { 
                        GelView(fragments: results) 
                    }
                    .padding()
                    .background(Color(.systemBackground))
                    .cornerRadius(8)
                    
                    NavigationLink("Ligation Analysis") { 
                        LigationCompatibilityView(fragments: results) 
                    }
                    .padding()
                    .background(Color(.systemBackground))
                    .cornerRadius(8)
                }
                .padding(.horizontal)
            }
        }
        .padding(.vertical)
    }
    
    private var selectedEnzymesView: some View {
        let maxDisplay = 6
        let sortedSelected = Array(selected).sorted { $0.name < $1.name }
        let displayedEnzymes = sortedSelected.prefix(maxDisplay)
        let remainingCount = selected.count - maxDisplay
        
        return VStack(alignment: .leading, spacing: 8) {
            Divider()
            
            Text("Selected Enzymes:")
                .font(.caption)
                .foregroundColor(.genomancerSecondaryText)
            
            FlowLayout(spacing: 8) {
                ForEach(Array(displayedEnzymes), id: \.self) { enzyme in
                    EnzymeChip(enzyme: enzyme, onRemove: {
                        selected.remove(enzyme)
                    })
                }
                
                if remainingCount > 0 {
                    Text("+\(remainingCount) more")
                        .font(.caption)
                        .padding(.horizontal, 10)
                        .padding(.vertical, 6)
                        .background(Color.genomancerSecondaryBackground)
                        .foregroundColor(.genomancerSecondaryText)
                        .cornerRadius(12)
                }
            }
        }
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
        // Enable sequence retrieval so fragments contain their DNA sequences
        results = engine.digest(options: .init(circular: circular, returnSequences: true))
    }

    func parseFASTA(_ s: String) -> String {
        s.split(separator: "\n").filter{ !$0.hasPrefix(">") }
            .joined().replacingOccurrences(of: "\\s", with: "", options: .regularExpression)
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

// MARK: - Enzyme Chip Component
struct EnzymeChip: View {
    let enzyme: Enzyme
    let onRemove: () -> Void
    
    var body: some View {
        HStack(spacing: 4) {
            Text(enzyme.name)
                .font(.caption)
                .lineLimit(1)
            
            Button(action: onRemove) {
                Image(systemName: "xmark.circle.fill")
                    .font(.system(size: 14))
                    .foregroundColor(.white.opacity(0.8))
            }
            .buttonStyle(.plain)
        }
        .padding(.horizontal, 10)
        .padding(.vertical, 6)
        .background(Color.genomancerSecondary)
        .foregroundColor(.genomancerText)
        .cornerRadius(12)
    }
}

// MARK: - Flow Layout
struct FlowLayout: Layout {
    var spacing: CGFloat = 8
    
    func sizeThatFits(proposal: ProposedViewSize, subviews: Subviews, cache: inout ()) -> CGSize {
        let result = FlowResult(
            in: proposal.replacingUnspecifiedDimensions().width,
            subviews: subviews,
            spacing: spacing
        )
        return result.size
    }
    
    func placeSubviews(in bounds: CGRect, proposal: ProposedViewSize, subviews: Subviews, cache: inout ()) {
        let result = FlowResult(
            in: bounds.width,
            subviews: subviews,
            spacing: spacing
        )
        for (index, subview) in subviews.enumerated() {
            subview.place(at: CGPoint(x: bounds.minX + result.positions[index].x, y: bounds.minY + result.positions[index].y), proposal: .unspecified)
        }
    }
    
    struct FlowResult {
        var size: CGSize = .zero
        var positions: [CGPoint] = []
        
        init(in maxWidth: CGFloat, subviews: Subviews, spacing: CGFloat) {
            var currentX: CGFloat = 0
            var currentY: CGFloat = 0
            var lineHeight: CGFloat = 0
            
            for subview in subviews {
                let subviewSize = subview.sizeThatFits(.unspecified)
                
                if currentX + subviewSize.width > maxWidth && currentX > 0 {
                    // Move to next line
                    currentX = 0
                    currentY += lineHeight + spacing
                    lineHeight = 0
                }
                
                positions.append(CGPoint(x: currentX, y: currentY))
                lineHeight = max(lineHeight, subviewSize.height)
                currentX += subviewSize.width + spacing
                size.width = max(size.width, currentX - spacing)
            }
            
            size.height = currentY + lineHeight
        }
    }
}


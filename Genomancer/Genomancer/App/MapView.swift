import SwiftUI
import WebKit
import DigestCore

struct MapView: View {
    let sequence: String
    let enzymes: [Enzyme]
    let circular: Bool
    
    @State private var svgContent: String = ""
    @State private var showCopiedAlert = false
    
    var body: some View {
        VStack {
            if svgContent.isEmpty {
                ProgressView("Generating map...")
                    .foregroundColor(.genomancerText)
                    .onAppear(perform: generateMap)
            } else {
                SVGWebView(svgContent: svgContent)
                    .navigationTitle("Plasmid Map")
                    .navigationBarTitleDisplayMode(.inline)
                    .accessibilityLabel("Plasmid map visualization")
                    .accessibilityHint(mapAccessibilityDescription)
                    .toolbar {
                        ToolbarItem(placement: .navigationBarTrailing) {
                            Button(action: copySVG) {
                                Label("Copy SVG", systemImage: "doc.on.doc")
                            }
                            .accessibilityLabel("Copy SVG to clipboard")
                        }
                    }
                    .alert("SVG Copied", isPresented: $showCopiedAlert) {
                        Button("OK", role: .cancel) { }
                    } message: {
                        Text("SVG code copied to clipboard")
                    }
            }
        }
        .background(Color.genomancerBackground)
    }
    
    private func copySVG() {
        UIPasteboard.general.string = svgContent
        showCopiedAlert = true
    }
    
    private func generateMap() {
        let engine = DigestEngine(sequence: sequence, enzymes: enzymes)
        let sites = engine.cutSites()
        svgContent = plasmidSVG(dnaLength: sequence.count, sites: sites)
    }
    
    private var mapAccessibilityDescription: String {
        let engine = DigestEngine(sequence: sequence, enzymes: enzymes)
        let sites = engine.cutSites()
        let enzymeNames = Set(sites.map { $0.enzyme.name })
        return "Plasmid map for \(circular ? "circular" : "linear") DNA, \(sequence.count) base pairs, showing \(sites.count) cut sites from \(enzymeNames.count) enzymes: \(enzymeNames.sorted().joined(separator: ", "))"
    }
}

/// WebView wrapper to display SVG content
struct SVGWebView: UIViewRepresentable {
    let svgContent: String
    
    func makeUIView(context: Context) -> WKWebView {
        let webView = WKWebView()
        webView.isOpaque = false
        webView.backgroundColor = .clear
        return webView
    }
    
    func updateUIView(_ webView: WKWebView, context: Context) {
        let html = """
        <!DOCTYPE html>
        <html>
        <head>
            <meta name="viewport" content="width=device-width, initial-scale=1.0, maximum-scale=1.0, user-scalable=no">
            <style>
                body {
                    margin: 0;
                    padding: 0;
                    display: flex;
                    justify-content: center;
                    align-items: center;
                    min-height: 100vh;
                }
                svg {
                    max-width: 100%;
                    height: auto;
                }
            </style>
        </head>
        <body>
            \(svgContent)
        </body>
        </html>
        """
        
        webView.loadHTMLString(html, baseURL: nil)
    }
}


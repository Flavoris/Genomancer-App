import SwiftUI
import WebKit
import DigestCore

struct MapView: View {
    let sequence: String
    let enzymes: [Enzyme]
    let circular: Bool
    
    @State private var svgContent: String = ""
    
    var body: some View {
        VStack {
            if svgContent.isEmpty {
                ProgressView("Generating map...")
                    .onAppear(perform: generateMap)
            } else {
                SVGWebView(svgContent: svgContent)
                    .navigationTitle("Plasmid Map")
                    .navigationBarTitleDisplayMode(.inline)
            }
        }
    }
    
    private func generateMap() {
        let engine = DigestEngine(sequence: sequence, enzymes: enzymes)
        let sites = engine.cutSites()
        svgContent = plasmidSVG(dnaLength: sequence.count, sites: sites)
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


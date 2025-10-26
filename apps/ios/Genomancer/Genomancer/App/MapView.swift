import SwiftUI
import WebKit
import DigestCore

struct MapView: View {
    let sequence: String
    let enzymes: [Enzyme]
    let circular: Bool
    
    @State private var svgContent: String = ""
    @State private var cutSites: [CutSite] = []
    
    var body: some View {
        VStack {
            if svgContent.isEmpty {
                ProgressView("Generating map...")
                    .foregroundColor(.genomancerText)
                    .onAppear(perform: generateMap)
            } else {
                InteractivePlasmidMapView(
                    svgContent: svgContent,
                    cutSites: cutSites,
                    sequence: sequence,
                    dnaLength: sequence.count
                )
                .navigationTitle("Plasmid Map")
                .navigationBarTitleDisplayMode(.inline)
                .toolbarColorScheme(.dark, for: .navigationBar)
                .accessibilityLabel("Plasmid map visualization")
                .accessibilityHint(mapAccessibilityDescription)
            }
        }
        .background(Color.genomancerBackground)
    }
    
    private func generateMap() {
        let engine = DigestEngine(sequence: sequence, enzymes: enzymes)
        let sites = engine.cutSites()
        cutSites = sites
        svgContent = plasmidSVG(dnaLength: sequence.count, sites: sites)
    }
    
    private var mapAccessibilityDescription: String {
        let engine = DigestEngine(sequence: sequence, enzymes: enzymes)
        let sites = engine.cutSites()
        let enzymeNames = Set(sites.map { $0.enzyme.name })
        return "Plasmid map for \(circular ? "circular" : "linear") DNA, \(sequence.count) base pairs, showing \(sites.count) cut sites from \(enzymeNames.count) enzymes: \(enzymeNames.sorted().joined(separator: ", "))"
    }
}

/// Interactive plasmid map with tappable cut sites
struct InteractivePlasmidMapView: View {
    let svgContent: String
    let cutSites: [CutSite]
    let sequence: String
    let dnaLength: Int
    
    @State private var showHint = true
    @State private var selectedCluster: ClusterWrapper?
    @State private var selectedSingleSite: IndexedCutSite?
    
    var body: some View {
        GeometryReader { geometry in
            let clusters = groupCutSitesByProximity(in: geometry.size)
            
            ZStack {
                // Background SVG
                SVGWebView(svgContent: svgContent)
                    .onAppear {
                        // Auto-dismiss hint after 5 seconds
                        DispatchQueue.main.asyncAfter(deadline: .now() + 5) {
                            withAnimation {
                                showHint = false
                            }
                        }
                    }
                
                // Overlay interactive buttons at cut site positions
                ForEach(clusters.indices, id: \.self) { clusterIndex in
                    let cluster = clusters[clusterIndex]
                    let position = calculateClusterCenter(cluster: cluster, in: geometry.size)
                    let isCluster = cluster.count > 1
                    
                    Button(action: {
                        if isCluster {
                            selectedCluster = ClusterWrapper(sites: cluster)
                        } else {
                            // Single site - show as sheet
                            selectedSingleSite = cluster[0]
                        }
                    }) {
                        ZStack {
                            // Pulsing outer ring for better visibility
                            Circle()
                                .stroke(Color.genomancerAccent.opacity(0.3), lineWidth: 2)
                                .frame(width: 44, height: 44)
                                .scaleEffect(showHint ? 1.2 : 1.0)
                                .animation(.easeInOut(duration: 1.5).repeatForever(autoreverses: true), value: showHint)
                            
                            // Inner filled circle with glow effect
                            Circle()
                                .fill(Color.genomancerAccent.opacity(0.5))
                                .frame(width: 32, height: 32)
                                .shadow(color: Color.genomancerAccent.opacity(0.6), radius: 8, x: 0, y: 0)
                            
                            // Icon to indicate tappability
                            if isCluster {
                                // Show count badge for clusters
                                ZStack {
                                    Image(systemName: "info.circle.fill")
                                        .font(.system(size: 14))
                                        .foregroundColor(.white)
                                        .shadow(color: .black.opacity(0.3), radius: 2, x: 0, y: 1)
                                    
                                    // Badge with count
                                    Text("\(cluster.count)")
                                        .font(.system(size: 10, weight: .bold))
                                        .foregroundColor(.white)
                                        .padding(4)
                                        .background(Circle().fill(Color.red))
                                        .offset(x: 12, y: -12)
                                }
                            } else {
                                Image(systemName: "info.circle.fill")
                                    .font(.system(size: 14))
                                    .foregroundColor(.white)
                                    .shadow(color: .black.opacity(0.3), radius: 2, x: 0, y: 1)
                            }
                        }
                    }
                    .buttonStyle(PlainButtonStyle())
                    .position(position)
                }
                
                // Hint overlay at the bottom
                if showHint && !cutSites.isEmpty {
                    VStack {
                        Spacer()
                        HStack {
                            Image(systemName: "hand.tap.fill")
                                .foregroundColor(.genomancerAccent)
                            Text("Tap any cut site to view details")
                                .font(.caption)
                                .foregroundColor(.genomancerText)
                            Button(action: {
                                withAnimation {
                                    showHint = false
                                }
                            }) {
                                Image(systemName: "xmark.circle.fill")
                                    .foregroundColor(.genomancerSecondaryText)
                            }
                        }
                        .padding()
                        .background(
                            RoundedRectangle(cornerRadius: 12)
                                .fill(Color.genomancerCardBackground.opacity(0.95))
                        )
                        .padding()
                    }
                    .transition(.move(edge: .bottom).combined(with: .opacity))
                }
            }
            .sheet(item: $selectedCluster) { clusterWrapper in
                ClusterSelectionView(cluster: clusterWrapper.sites, sequence: sequence)
                    .presentationDetents([.medium, .large])
                    .presentationDragIndicator(.visible)
            }
            .sheet(item: $selectedSingleSite) { indexedSite in
                NavigationStack {
                    CutSiteDetailView(
                        cutSite: indexedSite.site,
                        index: indexedSite.index,
                        sequence: sequence
                    )
                    .toolbar {
                        ToolbarItem(placement: .navigationBarTrailing) {
                            Button("Done") {
                                selectedSingleSite = nil
                            }
                        }
                    }
                }
                .presentationDetents([.medium, .large])
                .presentationDragIndicator(.visible)
            }
        }
    }
    
    // Helper struct to track cut sites with their original indices
    struct IndexedCutSite: Identifiable, Hashable {
        let id = UUID()
        let site: CutSite
        let index: Int
        let position: CGPoint
        
        // Implement Hashable based on the unique id
        func hash(into hasher: inout Hasher) {
            hasher.combine(id)
        }
        
        static func == (lhs: IndexedCutSite, rhs: IndexedCutSite) -> Bool {
            lhs.id == rhs.id
        }
    }
    
    // Helper struct to wrap a cluster for sheet presentation
    struct ClusterWrapper: Identifiable {
        let id = UUID()
        let sites: [IndexedCutSite]
    }
    
    /// Group cut sites that are close together (within 35 points)
    private func groupCutSitesByProximity(in size: CGSize) -> [[IndexedCutSite]] {
        let threshold: CGFloat = 35.0
        let indexedSites = cutSites.enumerated().map { index, site in
            IndexedCutSite(
                site: site,
                index: index,
                position: calculateCutSitePosition(for: site, in: size)
            )
        }
        
        var clusters: [[IndexedCutSite]] = []
        var used = Set<Int>()
        
        for i in 0..<indexedSites.count {
            if used.contains(i) { continue }
            
            var cluster = [indexedSites[i]]
            used.insert(i)
            
            for j in (i+1)..<indexedSites.count {
                if used.contains(j) { continue }
                
                let distance = sqrt(
                    pow(indexedSites[i].position.x - indexedSites[j].position.x, 2) +
                    pow(indexedSites[i].position.y - indexedSites[j].position.y, 2)
                )
                
                if distance < threshold {
                    cluster.append(indexedSites[j])
                    used.insert(j)
                }
            }
            
            clusters.append(cluster)
        }
        
        return clusters
    }
    
    /// Calculate the center position for a cluster of cut sites
    private func calculateClusterCenter(cluster: [IndexedCutSite], in size: CGSize) -> CGPoint {
        let avgX = cluster.map { $0.position.x }.reduce(0, +) / CGFloat(cluster.count)
        let avgY = cluster.map { $0.position.y }.reduce(0, +) / CGFloat(cluster.count)
        return CGPoint(x: avgX, y: avgY)
    }
    
    /// Calculate the screen position for a cut site marker
    private func calculateCutSitePosition(for site: CutSite, in size: CGSize) -> CGPoint {
        // SVG viewBox is "-50 -50 700 700"
        let viewBoxMinX: CGFloat = -50
        let viewBoxMinY: CGFloat = -50
        let viewBoxWidth: CGFloat = 700
        let viewBoxHeight: CGFloat = 700
        
        // Center of the plasmid circle in SVG coordinates
        let centerX: CGFloat = 300
        let centerY: CGFloat = 300
        let radius: CGFloat = 220
        
        // Calculate angle for this cut site (same as in Render.swift)
        let angle = (Double(site.position) / Double(dnaLength)) * 2.0 * .pi - .pi / 2.0
        
        // Calculate position on the circle in SVG coordinates
        let svgX = centerX + radius * CGFloat(cos(angle))
        let svgY = centerY + radius * CGFloat(sin(angle))
        
        // The SVG maintains aspect ratio and is centered in the view
        // Calculate the actual rendered size
        let aspectRatio = viewBoxWidth / viewBoxHeight
        let screenAspectRatio = size.width / size.height
        
        var renderedWidth: CGFloat
        var renderedHeight: CGFloat
        var offsetX: CGFloat
        var offsetY: CGFloat
        
        if screenAspectRatio > aspectRatio {
            // Screen is wider than SVG - SVG is constrained by height
            renderedHeight = size.height
            renderedWidth = renderedHeight * aspectRatio
            offsetX = (size.width - renderedWidth) / 2
            offsetY = 0
        } else {
            // Screen is taller than SVG - SVG is constrained by width
            renderedWidth = size.width
            renderedHeight = renderedWidth / aspectRatio
            offsetX = 0
            offsetY = (size.height - renderedHeight) / 2
        }
        
        // Scale from SVG coordinates to screen coordinates
        let scale = renderedWidth / viewBoxWidth
        
        // Convert SVG coordinates to screen coordinates
        let screenX = offsetX + (svgX - viewBoxMinX) * scale
        let screenY = offsetY + (svgY - viewBoxMinY) * scale
        
        return CGPoint(x: screenX, y: screenY)
    }
}

/// Sheet view for selecting a cut site from a cluster
struct ClusterSelectionView: View {
    typealias IndexedCutSite = InteractivePlasmidMapView.IndexedCutSite
    
    let cluster: [IndexedCutSite]
    let sequence: String
    
    @Environment(\.dismiss) private var dismiss
    
    var body: some View {
        NavigationStack {
            List {
                Section {
                    ForEach(cluster.sorted(by: { $0.site.position < $1.site.position })) { indexedSite in
                        NavigationLink(destination: CutSiteDetailView(
                            cutSite: indexedSite.site,
                            index: indexedSite.index,
                            sequence: sequence
                        )) {
                            HStack {
                                VStack(alignment: .leading, spacing: 4) {
                                    Text(indexedSite.site.enzyme.name)
                                        .font(.headline)
                                        .foregroundColor(.genomancerAccent)
                                    
                                    Text("\(indexedSite.site.position) bp")
                                        .font(.caption)
                                        .foregroundColor(.secondary)
                                    
                                    if let overhangType = indexedSite.site.enzyme.overhangType {
                                        Text(overhangTypeDescription(overhangType))
                                            .font(.caption2)
                                            .foregroundColor(.genomancerInfo)
                                    }
                                }
                                
                                Spacer()
                                
                                Image(systemName: "chevron.right")
                                    .font(.caption)
                                    .foregroundColor(.secondary)
                            }
                            .padding(.vertical, 4)
                        }
                    }
                } header: {
                    Text("Select Cut Site")
                } footer: {
                    Text("\(cluster.count) cut sites at this location")
                }
            }
            .navigationTitle("Clustered Cut Sites")
            .navigationBarTitleDisplayMode(.inline)
            .toolbar {
                ToolbarItem(placement: .navigationBarTrailing) {
                    Button("Done") {
                        dismiss()
                    }
                }
            }
        }
    }
    
    private func overhangTypeDescription(_ type: OverhangType) -> String {
        switch type {
        case .blunt:
            return "Blunt end"
        case .fivePrime:
            return "5' overhang"
        case .threePrime:
            return "3' overhang"
        case .unknown:
            return "Unknown"
        }
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


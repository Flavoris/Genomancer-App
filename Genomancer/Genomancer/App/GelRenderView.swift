import SwiftUI
import DigestCore

/// A minimal SwiftUI renderer for a single gel lane.
/// For now, it only draws a background and a centered lane with vertical padding.
struct GelRenderView: View {
    enum GelStyle { case photo, cartoon }

    let fragments: [Fragment]
    let gelPercent: Double
    let showTicks: Bool
    let style: GelStyle

    init(fragments: [Fragment], gelPercent: Double = 2.0, showTicks: Bool = true, style: GelStyle = .photo) {
        self.fragments = fragments
        self.gelPercent = gelPercent
        self.showTicks = showTicks
        self.style = style
    }

    var body: some View {
        GeometryReader { geometry in
            ZStack {
                Color.genomancerBackground
                Canvas { context, size in
                    // Layout constants and style tuning
                    let verticalPadding: CGFloat = 24
                    let laneWidth: CGFloat = 110

                    let laneHeight = max(0, size.height - (verticalPadding * 2))
                    // Center the lane in the view (accounting for label space on left)
                    let labelSpace: CGFloat = 80
                    let laneX = labelSpace + (size.width - labelSpace - laneWidth) / 2
                    let laneY = verticalPadding
                    let laneRect = CGRect(x: laneX, y: laneY, width: laneWidth, height: laneHeight)

                    // Lane background (subtle), with a thin border
                    let laneBackground = Color.white.opacity(0.06)
                    let laneBorder = Color.white.opacity(0.08)

                    let path = Path(roundedRect: laneRect, cornerRadius: 6)
                    context.fill(path, with: .color(laneBackground))
                    // Draw fragment bands using semi-log size→y mapping
                    drawBands(context: context, size: size, laneRect: laneRect, gelPercent: gelPercent)
                    context.stroke(path, with: .color(laneBorder), lineWidth: 1)
                }
                .compositingGroup()
                .drawingGroup()
                
                // Overlay text labels outside Canvas
                if showTicks {
                    tickLabelsOverlay(geometry: geometry)
                }
            }
        }
        .accessibilityLabel("Gel electrophoresis lane")
        .accessibilityValue(laneAccessibilityDescription)
    }

    private var laneAccessibilityDescription: String {
        let lengths = fragments.map { $0.length }.sorted(by: >)
        if lengths.isEmpty { return "No fragments" }
        let top = lengths.prefix(5).map { "\($0) base pairs" }.joined(separator: ", ")
        let suffix = lengths.count > 5 ? " and \(lengths.count - 5) more" : ""
        return "Showing \(lengths.count) fragments: \(top)\(suffix)"
    }

    // MARK: - Rendering helpers
    
    /// Creates text label overlays for fragment bands
    @ViewBuilder
    private func tickLabelsOverlay(geometry: GeometryProxy) -> some View {
        let lengths = fragments.map { $0.length }
        if let minLength = lengths.min(), let maxLength = lengths.max(), minLength > 0, maxLength > 0 {
            let size = geometry.size
            let verticalPadding: CGFloat = 24
            let laneWidth: CGFloat = 110
            
            let laneHeight = max(0, size.height - (verticalPadding * 2))
            // Center the lane in the view (accounting for label space on left)
            let labelSpace: CGFloat = 80
            let laneX = labelSpace + (size.width - labelSpace - laneWidth) / 2
            let laneY = verticalPadding
            let laneRect = CGRect(x: laneX, y: laneY, width: laneWidth, height: laneHeight)
            
            let maxCoreHeight: CGFloat = 8
            let mappingTop = laneRect.minY + maxCoreHeight / 2
            let mappingHeight = max(0, laneRect.height - maxCoreHeight)
            
            let labelX = laneRect.minX - 20
            let labelColor = Color.white.opacity(0.6)
            
            // Get unique fragment sizes, sorted descending
            let uniqueLengths = Array(Set(lengths)).sorted(by: >)
            
            ForEach(uniqueLengths, id: \.self) { bp in
                let y = GelLayout.y(bp: bp, minBP: minLength, maxBP: maxLength, top: mappingTop, height: mappingHeight, gelPercent: gelPercent)
                if y >= laneRect.minY && y <= laneRect.maxY {
                    Text("\(bp) bp")
                        .font(.system(size: 13, weight: .medium, design: .rounded))
                        .foregroundStyle(labelColor)
                        .position(x: labelX, y: y)
                }
            }
        }
    }

    /// Draws fragment bands within the provided lane rectangle using a semi-log mapping.
    private func drawBands(context: GraphicsContext, size: CGSize, laneRect: CGRect, gelPercent: Double) {
        guard !fragments.isEmpty else { return }

        // Compute min and max base-pair sizes from fragments
        let lengths = fragments.map { $0.length }
        guard let minLength = lengths.min(), let maxLength = lengths.max(), minLength > 0, maxLength > 0 else { return }

        // Inner padding within lane for band width
        let horizontalPadding: CGFloat = 8
        let innerX = laneRect.minX + horizontalPadding
        let innerWidth = max(0, laneRect.width - (horizontalPadding * 2))

        // Mapping region for y-centers: keep room for max band thickness at top/bottom
        let maxCoreHeight: CGFloat = 8
        let mappingTop = laneRect.minY + maxCoreHeight / 2
        let mappingHeight = max(0, laneRect.height - maxCoreHeight)

        // Sort descending by length so larger fragments (higher on gel) draw first
        let sortedFragments = fragments.sorted { $0.length > $1.length }

        // Offscreen layer for optional subtle blur over the entire lane bands
        let laneClipPath = Path(roundedRect: laneRect, cornerRadius: 6)
        context.drawLayer { layer in
            // Keep drawing confined to the lane
            layer.clip(to: laneClipPath)

            // Small Gaussian blur to gently soften the overall lane content
            // Keep radius tiny to avoid heavy performance impact
            layer.addFilter(.blur(radius: 0.9))

            for fragment in sortedFragments {
                let bp = max(1, fragment.length)

                // y-center using semi-log mapping with gel percentage
                let yCenter = GelLayout.y(bp: bp, minBP: minLength, maxBP: maxLength, top: mappingTop, height: mappingHeight, gelPercent: gelPercent)

                // Band thickness heuristic
                var coreHeight = max(4.0, 8.0 - CGFloat(sqrt(Double(bp))) / 30.0)
                var haloOpacity: CGFloat = 0.35
                var haloRadius: CGFloat = 10
                var coreShadowOpacity: CGFloat = 0.6

                switch style {
                case .photo:
                    haloOpacity = 0.40
                    haloRadius = 12
                    coreShadowOpacity = 0.7
                case .cartoon:
                    coreHeight += 1.0
                    haloOpacity = 0.25
                    haloRadius = 7
                    coreShadowOpacity = 0.55
                }
                let cornerRadius = coreHeight / 2.0

                // Center the band vertically on yCenter
                let bandRect = CGRect(
                    x: innerX,
                    y: yCenter - coreHeight / 2.0,
                    width: innerWidth,
                    height: coreHeight
                )
                let bandPath = Path(roundedRect: bandRect, cornerRadius: cornerRadius)

                // Pass 1: Halo
                layer.drawLayer { halo in
                    halo.addFilter(.shadow(color: .white.opacity(haloOpacity), radius: haloRadius, x: 0, y: 0))
                    halo.fill(bandPath, with: .color(.white.opacity(haloOpacity)))
                }

                // Pass 2: Core
                layer.drawLayer { core in
                    core.addFilter(.shadow(color: .white.opacity(coreShadowOpacity), radius: 2, x: 0, y: 0))
                    core.fill(bandPath, with: .color(.white))
                }
            }
        }
    }
}


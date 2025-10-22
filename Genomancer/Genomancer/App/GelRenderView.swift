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
    
    /// Computes which label indices should be placed on the right side to avoid overlap
    private func computeLabelsOnRight(labelData: [(bp: Int, y: CGFloat)], minSpacing: CGFloat) -> Set<Int> {
        var labelsOnRight = Set<Int>()
        
        for i in 0..<labelData.count {
            for j in (i+1)..<labelData.count {
                let spacing = abs(labelData[i].y - labelData[j].y)
                if spacing < minSpacing {
                    // Move the lower (smaller fragment) to the right
                    labelsOnRight.insert(j)
                }
            }
        }
        
        return labelsOnRight
    }
    
    /// Creates text label overlays for fragment bands
    @ViewBuilder
    private func tickLabelsOverlay(geometry: GeometryProxy) -> some View {
        let lengths = fragments.map { $0.length }
        if let actualMinBP = lengths.min(), let actualMaxBP = lengths.max(), actualMinBP > 0, actualMaxBP > 0 {
            // Use same dynamic range calculation as drawBands for consistent positioning
            let logMin = log10(Double(actualMinBP))
            let logMax = log10(Double(actualMaxBP))
            let logRange = logMax - logMin
            
            // Same padding factor as drawBands
            let paddingFactor = 0.25
            let expandedLogMin = logMin - (logRange * paddingFactor)
            let expandedLogMax = logMax + (logRange * paddingFactor)
            
            let theoreticalMinBP = max(1, Int(pow(10.0, expandedLogMin)))
            let theoreticalMaxBP = Int(pow(10.0, expandedLogMax))
            
            let size = geometry.size
            let verticalPadding: CGFloat = 24
            let laneWidth: CGFloat = 110
            
            let laneHeight = max(0, size.height - (verticalPadding * 2))
            // Center the lane in the view (accounting for label space on left)
            let labelSpace: CGFloat = 80
            let laneX = labelSpace + (size.width - labelSpace - laneWidth) / 2
            let laneY = verticalPadding
            let laneRect = CGRect(x: laneX, y: laneY, width: laneWidth, height: laneHeight)
            
            // Use same mapping as drawBands for consistent positioning
            let maxCoreHeight: CGFloat = 8
            let mappingTop = laneRect.minY + maxCoreHeight / 2
            let mappingHeight = max(0, laneRect.height - maxCoreHeight)
            
            let leftLabelX = laneRect.minX - 20
            let rightLabelX = laneRect.maxX + 20
            let labelColor = Color.white.opacity(0.6)
            
            // Get unique fragment sizes, sorted descending with their y positions
            let uniqueLengths = Array(Set(lengths)).sorted(by: >)
            let labelData: [(bp: Int, y: CGFloat)] = uniqueLengths.compactMap { bp in
                let y = GelLayout.y(bp: bp, minBP: theoreticalMinBP, maxBP: theoreticalMaxBP, top: mappingTop, height: mappingHeight, gelPercent: gelPercent)
                if y >= laneRect.minY && y <= laneRect.maxY {
                    return (bp, y)
                }
                return nil
            }
            
            // Determine which labels should go on the right to avoid overlap
            let minLabelSpacing: CGFloat = 18  // Minimum pixels between label centers
            let labelsOnRight = computeLabelsOnRight(labelData: labelData, minSpacing: minLabelSpacing)
            
            ForEach(Array(labelData.enumerated()), id: \.element.bp) { index, data in
                let isOnRight = labelsOnRight.contains(index)
                let xPos = isOnRight ? rightLabelX : leftLabelX
                
                Text("\(data.bp) bp")
                    .font(.system(size: 13, weight: .medium, design: .rounded))
                    .foregroundStyle(labelColor)
                    .position(x: xPos, y: data.y)
            }
        }
    }

    /// Draws fragment bands within the provided lane rectangle using a semi-log mapping.
    private func drawBands(context: GraphicsContext, size: CGSize, laneRect: CGRect, gelPercent: Double) {
        guard !fragments.isEmpty else { return }

        // Compute actual fragment range, then expand it to provide margin at both ends
        let lengths = fragments.map { $0.length }
        guard let actualMinBP = lengths.min(), let actualMaxBP = lengths.max(), actualMinBP > 0, actualMaxBP > 0 else { return }
        
        // Expand range in log space to leave space at top and bottom of gel
        // This ensures fragments move when gel% changes while keeping appropriate spacing
        let logMin = log10(Double(actualMinBP))
        let logMax = log10(Double(actualMaxBP))
        let logRange = logMax - logMin
        
        // Add 25% padding on each end in log space for realistic gel margins
        let paddingFactor = 0.25
        let expandedLogMin = logMin - (logRange * paddingFactor)
        let expandedLogMax = logMax + (logRange * paddingFactor)
        
        let theoreticalMinBP = max(1, Int(pow(10.0, expandedLogMin)))
        let theoreticalMaxBP = Int(pow(10.0, expandedLogMax))

        // Inner padding within lane for band width
        let horizontalPadding: CGFloat = 8
        let innerX = laneRect.minX + horizontalPadding
        let innerWidth = max(0, laneRect.width - (horizontalPadding * 2))

        // Mapping region: use full lane height since expanded range provides the margins
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

                // y-center using semi-log mapping with gel percentage against theoretical range
                let yCenter = GelLayout.y(bp: bp, minBP: theoreticalMinBP, maxBP: theoreticalMaxBP, top: mappingTop, height: mappingHeight, gelPercent: gelPercent)

                // Band thickness heuristic with gel% adjustment
                // Minor polish: make thin bands at higher gel% 
                let gelThicknessFactor = max(0.8, 1.1 - 0.08 * (gelPercent - 1.5))
                var coreHeight = max(4.0, 8.0 - CGFloat(sqrt(Double(bp))) / 30.0) * gelThicknessFactor
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


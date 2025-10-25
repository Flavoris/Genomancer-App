import Foundation

public func plasmidSVG(dnaLength: Int, sites: [CutSite]) -> String {
    let centerX = 300.0
    let centerY = 300.0
    let radius = 220.0
    let markerLength = 25.0
    
    var cutSiteMarkers = ""
    
    // Generate markers and labels for each cut site
    for site in sites {
        let angle = (Double(site.position) / Double(dnaLength)) * 2.0 * .pi - .pi / 2.0
        
        // Calculate positions for the marker line (from circle to outer point)
        let innerX = centerX + radius * cos(angle)
        let innerY = centerY + radius * sin(angle)
        let outerX = centerX + (radius + markerLength) * cos(angle)
        let outerY = centerY + (radius + markerLength) * sin(angle)
        
        // Draw marker line
        cutSiteMarkers += """
          <line x1="\(innerX)" y1="\(innerY)" x2="\(outerX)" y2="\(outerY)" stroke="#F28C40" stroke-width="2"/>
          <circle cx="\(innerX)" cy="\(innerY)" r="4" fill="#F28C40"/>
        
        """
        
        // Calculate position for text label (further out from marker)
        let labelDistance = radius + markerLength + 35.0
        let labelX = centerX + labelDistance * cos(angle)
        let labelY = centerY + labelDistance * sin(angle)
        
        // Adjust text anchor based on position around circle
        let textAnchor = labelX > centerX ? "start" : (labelX < centerX ? "end" : "middle")
        
        // Add enzyme name and position labels
        cutSiteMarkers += """
          <text x="\(labelX)" y="\(labelY - 10)" font-family="monospace" font-size="16" font-weight="bold" text-anchor="\(textAnchor)" fill="#F28C40">
            \(site.enzyme.name)
          </text>
          <text x="\(labelX)" y="\(labelY + 12)" font-family="monospace" font-size="15" text-anchor="\(textAnchor)" fill="#FFFFFF">
            \(site.position) bp
          </text>
        
        """
    }
    
    return """
    <svg viewBox="-50 -50 700 700" xmlns="http://www.w3.org/2000/svg">
      <circle cx="300" cy="300" r="220" fill="none" stroke="#59A673" stroke-width="3"/>
      <text x="300" y="305" font-family="monospace" font-size="16" font-weight="bold" text-anchor="middle" fill="#FFFFFF">
        \(dnaLength) bp
      </text>
      \(cutSiteMarkers)
    </svg>
    """
}


import Foundation

public func plasmidSVG(dnaLength: Int, sites: [CutSite]) -> String {
    """
    <svg viewBox="0 0 600 600" xmlns="http://www.w3.org/2000/svg">
      <circle cx="300" cy="300" r="220" fill="none" stroke="black" stroke-width="2"/>
      <text x="300" y="30" font-family="monospace" font-size="14" text-anchor="middle">
        \(dnaLength) bp
      </text>
    </svg>
    """
}


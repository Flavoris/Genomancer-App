import Foundation
import CoreGraphics

/// Pure layout utilities for mapping DNA fragment sizes to y-coordinates.
/// Kept free of SwiftUI for easy unit testing.
public struct GelLayout {
    /// Maps base pairs to a y position using a semi-log transform, adjusted by gel percentage.
    /// y(bp) = top + (1 - (log10(bp)-log10(minBP)) / (log10(maxBP)-log10(minBP))) * height
    /// - Parameters:
    ///   - bp: Fragment length in base pairs. Clamped to at least 1.
    ///   - minBP: Minimum bp in the lane. Clamped to at least 1.
    ///   - maxBP: Maximum bp in the lane. At least minBP+1 to avoid division by zero.
    ///   - top: The y position representing the largest fragment.
    ///   - height: The drawable height from top to bottom.
    ///   - gelPercent: Agarose gel percentage (0.8-3.0%). Lower % = better large fragment separation. Higher % = better small fragment separation.
    /// - Returns: A y coordinate where larger bp values are closer to `top`.
    public static func y(bp: Int, minBP: Int, maxBP: Int, top: CGFloat, height: CGFloat, gelPercent: Double = 2.0) -> CGFloat {
        let clampedBP = max(1, bp)
        let clampedMin = max(1, minBP)
        let clampedMax = max(clampedMin + 1, maxBP)

        let logMin = log10(Double(clampedMin))
        let logMax = log10(Double(clampedMax))
        let logBP = log10(Double(clampedBP))

        let denom = logMax - logMin
        if denom <= 0 { return top + height * 0.5 }

        // Base normalized position (0 = top/largest, 1 = bottom/smallest)
        let t = (logBP - logMin) / denom
        
        // Apply gel percentage effect:
        // Lower % (e.g., 0.8%) compresses small fragments, expands large fragments
        // Higher % (e.g., 3.0%) expands small fragments, compresses large fragments
        // We use a power function to simulate this: higher gel % increases the exponent
        let gelFactor = 0.5 + (gelPercent / 4.0) // Maps 0.8-3.0% to ~0.7-1.25
        let adjustedT = pow(t, gelFactor)
        
        let inverted = 1.0 - adjustedT
        return top + CGFloat(inverted) * height
    }
}



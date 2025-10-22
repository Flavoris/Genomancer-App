import Foundation
import CoreGraphics

/// Pure layout utilities for mapping DNA fragment sizes to y-coordinates.
/// Kept free of SwiftUI for easy unit testing.
public struct GelLayout {
    /// Maps base pairs to a y position using a semi-log transform.
    /// y(bp) = top + (1 - (log10(bp)-log10(minBP)) / (log10(maxBP)-log10(minBP))) * height
    /// - Parameters:
    ///   - bp: Fragment length in base pairs. Clamped to at least 1.
    ///   - minBP: Minimum bp in the lane. Clamped to at least 1.
    ///   - maxBP: Maximum bp in the lane. At least minBP+1 to avoid division by zero.
    ///   - top: The y position representing the largest fragment.
    ///   - height: The drawable height from top to bottom.
    /// - Returns: A y coordinate where larger bp values are closer to `top`.
    public static func y(bp: Int, minBP: Int, maxBP: Int, top: CGFloat, height: CGFloat) -> CGFloat {
        let clampedBP = max(1, bp)
        let clampedMin = max(1, minBP)
        let clampedMax = max(clampedMin + 1, maxBP)

        let logMin = log10(Double(clampedMin))
        let logMax = log10(Double(clampedMax))
        let logBP = log10(Double(clampedBP))

        let denom = logMax - logMin
        if denom <= 0 { return top + height * 0.5 }

        let t = (logBP - logMin) / denom
        let inverted = 1.0 - t
        return top + CGFloat(inverted) * height
    }
}



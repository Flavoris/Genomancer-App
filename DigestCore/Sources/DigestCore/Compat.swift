import Foundation

/// Determines if two fragment ends are compatible for ligation based on concrete sequence matching
/// - Parameters:
///   - a: First fragment end
///   - b: Second fragment end
///   - includeBlunt: Whether to consider blunt-blunt ends as compatible
///   - minOverhang: Minimum overhang length required for sticky ends
/// - Returns: True if the ends can ligate together
public func compatibleConcrete(a: EndInfo, b: EndInfo,
                               includeBlunt: Bool = false,
                               minOverhang: Int = 1) -> Bool {
    // sticky↔blunt always false
    let aSticky = (a.overhangLen > 0 && a.overhangSeq5to3?.isEmpty == false)
    let bSticky = (b.overhangLen > 0 && b.overhangSeq5to3?.isEmpty == false)
    if aSticky != bSticky { return false }
    if !aSticky && !bSticky { return includeBlunt } // blunt↔blunt
    
    if a.overhangType != b.overhangType { return false }
    if a.overhangLen < minOverhang || b.overhangLen < minOverhang { return false }
    if a.overhangLen != b.overhangLen { return false }
    
    guard let sa = a.overhangSeq5to3, let sb = b.overhangSeq5to3 else { return false }
    return sa.uppercased() == revcomp(sb).uppercased()
}


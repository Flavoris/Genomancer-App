import Foundation

public func asciiGel(fragments: [Int]) -> String {
    // naive scale mapping
    guard let maxbp = fragments.max() else { return "" }
    let lines = stride(from: maxbp, through: 0, by: -max(1, maxbp/20)).map { y -> String in
        let mark = fragments.contains { abs($0 - y) < max(1, maxbp/80) } ? "###" : " | "
        return String(format: "%6d %@", y, mark)
    }
    return lines.joined(separator: "\n")
}


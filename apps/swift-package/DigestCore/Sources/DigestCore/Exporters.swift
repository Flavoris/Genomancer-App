import Foundation

public func exportCSV(fragments: [Fragment]) -> String {
    var rows = ["start,end,length,left_overhang,right_overhang"]
    for f in fragments {
        let left = f.leftEnd.overhangSeq5to3 ?? ""
        let right = f.rightEnd.overhangSeq5to3 ?? ""
        rows.append("\(f.start),\(f.end),\(f.length),\(left),\(right)")
    }
    return rows.joined(separator: "\n")
}

public func exportGenBank(sequence: String, locus: String, features: [(Range<Int>, String)]) -> String {
    let locusName = locus.replacingOccurrences(of: "[^A-Za-z0-9_]", with: "_", options: .regularExpression)
    var out: [String] = []
    out.append(String(format:"LOCUS       %-16s %7d bp    DNA     PLN       01-JAN-2000", locusName, sequence.count))
    out.append("FEATURES             Location/Qualifiers")
    for feat in features {
        out.append(String(format:"     misc_feature    %d..%d", feat.0.lowerBound+1, feat.0.upperBound))
        out.append("                     /note=\"\(feat.1)\"")
    }
    out.append("ORIGIN")
    out.append(formatFASTAStyle(sequence))
    out.append("//")
    return out.joined(separator: "\n")
}

private func formatFASTAStyle(_ s: String) -> String {
    let arr = Array(s.lowercased())
    var lines: [String] = []
    var i = 0; var idx = 1
    while i < arr.count {
        let chunk = Array(arr[i..<min(i+60, arr.count)])
        let groups = stride(from: 0, to: chunk.count, by: 10).map { String(chunk[$0..<min($0+10, chunk.count)]) }
        lines.append(String(format:"%9d %@", idx, groups.joined(separator: " ")))
        i += 60; idx += 60
    }
    return lines.joined(separator: "\n")
}


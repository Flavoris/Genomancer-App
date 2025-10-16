import Foundation

public enum OverhangType: String, Codable {
    case fivePrime = "5'"
    case threePrime = "3'"
    case blunt
    case unknown
}

public struct Enzyme: Codable, Hashable {
    public let name: String
    public let site: String                // IUPAC motif
    public let cutIndexTop: Int?
    public let cutIndexBottom: Int?
    public let overhangType: OverhangType?
    public let notes: String?
    
    public init(
        name: String,
        site: String,
        cutIndexTop: Int?,
        cutIndexBottom: Int?,
        overhangType: OverhangType?,
        notes: String?
    ) {
        self.name = name
        self.site = site
        self.cutIndexTop = cutIndexTop
        self.cutIndexBottom = cutIndexBottom
        self.overhangType = overhangType
        self.notes = notes
    }
}


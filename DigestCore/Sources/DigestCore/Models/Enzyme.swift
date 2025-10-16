import Foundation

public enum OverhangType: String, Codable {
    case fivePrime = "5'"
    case threePrime = "3'"
    case blunt
    case unknown
}

public struct Enzyme: Codable, Hashable {
    public let name: String
    public let site: String                // IUPAC motif (stored uppercased)
    public let cutIndexTop: Int?
    public let cutIndexBottom: Int?
    public let overhangType: OverhangType?
    public let notes: String?
    
    // Performance optimization: cached motif as character array
    public let siteChars: [Character]
    
    public init(
        name: String,
        site: String,
        cutIndexTop: Int?,
        cutIndexBottom: Int?,
        overhangType: OverhangType?,
        notes: String?
    ) {
        self.name = name
        self.site = site.uppercased()  // Ensure motif is always uppercased
        self.cutIndexTop = cutIndexTop
        self.cutIndexBottom = cutIndexBottom
        self.overhangType = overhangType
        self.notes = notes
        self.siteChars = Array(self.site)
    }
    
    // Codable conformance with custom decoding
    enum CodingKeys: String, CodingKey {
        case name, site, cutIndexTop, cutIndexBottom, overhangType, notes
    }
    
    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        name = try container.decode(String.self, forKey: .name)
        let rawSite = try container.decode(String.self, forKey: .site)
        site = rawSite.uppercased()
        cutIndexTop = try container.decodeIfPresent(Int.self, forKey: .cutIndexTop)
        cutIndexBottom = try container.decodeIfPresent(Int.self, forKey: .cutIndexBottom)
        overhangType = try container.decodeIfPresent(OverhangType.self, forKey: .overhangType)
        notes = try container.decodeIfPresent(String.self, forKey: .notes)
        siteChars = Array(site)
    }
    
    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        try container.encode(name, forKey: .name)
        try container.encode(site, forKey: .site)
        try container.encodeIfPresent(cutIndexTop, forKey: .cutIndexTop)
        try container.encodeIfPresent(cutIndexBottom, forKey: .cutIndexBottom)
        try container.encodeIfPresent(overhangType, forKey: .overhangType)
        try container.encodeIfPresent(notes, forKey: .notes)
    }
    
    // Hashable conformance (exclude derived siteChars)
    public func hash(into hasher: inout Hasher) {
        hasher.combine(name)
        hasher.combine(site)
        hasher.combine(cutIndexTop)
        hasher.combine(cutIndexBottom)
        hasher.combine(overhangType)
        hasher.combine(notes)
    }
    
    public static func == (lhs: Enzyme, rhs: Enzyme) -> Bool {
        return lhs.name == rhs.name &&
               lhs.site == rhs.site &&
               lhs.cutIndexTop == rhs.cutIndexTop &&
               lhs.cutIndexBottom == rhs.cutIndexBottom &&
               lhs.overhangType == rhs.overhangType &&
               lhs.notes == rhs.notes
    }
}


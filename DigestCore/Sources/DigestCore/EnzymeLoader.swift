import Foundation

public enum EnzymeLoader {
    public static func load(from url: URL) throws -> [Enzyme] {
        let data = try Data(contentsOf: url)
        let dec = JSONDecoder()
        return try dec.decode([Enzyme].self, from: data)
    }
}


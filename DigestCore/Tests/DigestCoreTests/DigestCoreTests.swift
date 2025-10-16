import XCTest
@testable import DigestCore

final class DigestCoreTests: XCTestCase {
    
    func testDigestCoreInitialization() throws {
        let digestCore = DigestCore()
        XCTAssertNotNil(digestCore)
    }
    
    func testAnalyzeSequence() throws {
        let digestCore = DigestCore()
        let result = digestCore.analyzeSequence("ATCG")
        XCTAssertEqual(result, "Analyzing sequence: ATCG")
    }
}

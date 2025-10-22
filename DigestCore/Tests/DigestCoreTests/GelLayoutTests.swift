import XCTest
@testable import DigestCore

final class GelLayoutTests: XCTestCase {
    func testTopAtMaxBP() {
        let top: CGFloat = 10
        let height: CGFloat = 200
        let minBP = 100
        let maxBP = 5000
        let y = GelLayout.y(bp: maxBP, minBP: minBP, maxBP: maxBP, top: top, height: height)
        XCTAssertEqual(y, top, accuracy: 0.5)
    }

    func testBottomAtMinBP() {
        let top: CGFloat = 10
        let height: CGFloat = 200
        let minBP = 100
        let maxBP = 5000
        let y = GelLayout.y(bp: minBP, minBP: minBP, maxBP: maxBP, top: top, height: height)
        XCTAssertEqual(y, top + height, accuracy: 0.5)
    }

    func testMonotonicity() {
        let top: CGFloat = 0
        let height: CGFloat = 300
        let minBP = 100
        let maxBP = 10000
        let values = [120, 200, 500, 1000, 2000, 5000, 9000]
        let ys = values.map { GelLayout.y(bp: $0, minBP: minBP, maxBP: maxBP, top: top, height: height) }
        // Larger bp should map higher (smaller y)
        for i in 1..<ys.count {
            XCTAssertLessThanOrEqual(ys[i], ys[i-1] + 1e-6)
        }
    }
}



// swift-tools-version: 5.9
import PackageDescription

let package = Package(
    name: "DigestCore",
    platforms: [
        .iOS(.v16)
    ],
    products: [
        .library(
            name: "DigestCore",
            targets: ["DigestCore"]),
    ],
    targets: [
        .target(
            name: "DigestCore",
            dependencies: [],
            path: "Sources/DigestCore"),
        .testTarget(
            name: "DigestCoreTests",
            dependencies: ["DigestCore"],
            path: "Tests/DigestCoreTests"),
    ]
)

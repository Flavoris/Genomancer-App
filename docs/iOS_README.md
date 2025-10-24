# Genomancer iOS App

This is an iOS app for DNA digestion analysis, built with SwiftUI and a custom Swift package.

## Project Structure

```
Genomancer/
├── Genomancer.xcworkspace/          # Xcode workspace file
├── DigestCore/                      # Swift Package
│   ├── Package.swift
│   ├── Sources/DigestCore/
│   │   └── DigestCore.swift
│   └── Tests/DigestCoreTests/
│       └── DigestCoreTests.swift
└── Genomancer/                      # iOS App
    ├── Genomancer.xcodeproj/
    ├── App/
    │   ├── GenomancerApp.swift
    │   └── ContentView.swift
    └── Resources/
        ├── Assets.xcassets/
        └── Preview Content/
```

## Requirements

- Xcode 15.0+
- iOS 16.0+
- Swift 5.9+

## How to Build

1. Open `Genomancer.xcworkspace` in Xcode
2. Select the Genomancer scheme
3. Choose your target device or simulator
4. Press Cmd+R to build and run

## Features

- Basic SwiftUI interface
- Integration with DigestCore Swift package
- DNA sequence analysis functionality
- Clean, modular architecture

## Next Steps

The app currently has a basic interface. You can now:

1. Port your Python DNA analysis logic to the DigestCore Swift package
2. Enhance the UI with more sophisticated features
3. Add data persistence
4. Implement additional analysis tools

## Development Notes

- The app depends on the local DigestCore package
- All Swift files follow iOS development best practices
- The project is configured for iOS 16+ deployment

import SwiftUI

// Color theme matching the logo design - dark teal background with wizard hat and DNA helix
extension Color {
    // Primary brand colors - matching the logo's dark teal background and dark blue wizard hat
    static let genomancerPrimary = Color(red: 0.05, green: 0.15, blue: 0.20) // #0C2633 - exact logo background color
    static let genomancerSecondary = Color(red: 0.15, green: 0.30, blue: 0.45) // #264D73 - dark blue wizard hat
    static let genomancerAccent = Color(red: 0.95, green: 0.55, blue: 0.25) // #F28C40 - bright orange from DNA bases
    
    // DNA helix inspired colors - exact colors from the logo's DNA base pairs
    static let genomancerDNABase1 = Color(red: 0.95, green: 0.55, blue: 0.25) // #F28C40 - bright orange
    static let genomancerDNABase2 = Color(red: 0.95, green: 0.75, blue: 0.25) // #F2BF40 - golden yellow
    static let genomancerDNABase3 = Color(red: 0.35, green: 0.65, blue: 0.85) // #59A6D9 - medium blue
    static let genomancerDNABase4 = Color(red: 0.85, green: 0.45, blue: 0.35) // #D97359 - reddish-orange
    
    // DNA helix strand color - medium green from the logo
    static let genomancerDNAHelix = Color(red: 0.35, green: 0.65, blue: 0.45) // #59A673 - medium green helix
    
    // Background variations - creating depth with the dark teal theme
    static let genomancerBackground = genomancerPrimary
    static let genomancerCardBackground = Color(red: 0.08, green: 0.20, blue: 0.28) // Slightly lighter than primary
    static let genomancerSecondaryBackground = Color(red: 0.10, green: 0.25, blue: 0.32) // Even lighter
    
    // Surface colors for cards and overlays - matching the dark teal theme
    static let genomancerSurface = Color(red: 0.18, green: 0.35, blue: 0.45) // Light surface
    static let genomancerSurfaceVariant = Color(red: 0.22, green: 0.40, blue: 0.52) // Lighter variant
    
    // Text colors - high contrast for accessibility
    static let genomancerText = Color.white
    static let genomancerSecondaryText = Color(red: 0.85, green: 0.90, blue: 0.95) // Slightly off-white
    static let genomancerTertiaryText = Color(red: 0.70, green: 0.80, blue: 0.85) // Muted text
    
    // Status colors - using DNA base colors for consistency
    static let genomancerSuccess = genomancerDNAHelix // Success green (DNA helix color)
    static let genomancerWarning = genomancerDNABase2 // Warning golden yellow
    static let genomancerError = Color(red: 0.85, green: 0.35, blue: 0.35) // Error red
    static let genomancerInfo = genomancerDNABase3 // Info medium blue
}

// Custom button styles
struct GenomancerButtonStyle: ButtonStyle {
    func makeBody(configuration: Configuration) -> some View {
        configuration.label
            .foregroundColor(.genomancerText)
            .padding()
            .background(
                RoundedRectangle(cornerRadius: 12)
                    .fill(Color.genomancerSecondary)
                    .opacity(configuration.isPressed ? 0.8 : 1.0)
            )
            .scaleEffect(configuration.isPressed ? 0.95 : 1.0)
            .animation(.easeInOut(duration: 0.1), value: configuration.isPressed)
    }
}

struct GenomancerProminentButtonStyle: ButtonStyle {
    @Environment(\.isEnabled) private var isEnabled
    
    func makeBody(configuration: Configuration) -> some View {
        configuration.label
            .foregroundColor(isEnabled ? .white : Color.white.opacity(0.5))
            .padding()
            .background(
                RoundedRectangle(cornerRadius: 12)
                    .fill(isEnabled ? Color(red: 0.2, green: 0.65, blue: 0.3) : Color.gray.opacity(0.4))
                    .opacity(configuration.isPressed ? 0.8 : 1.0)
            )
            .scaleEffect(configuration.isPressed ? 0.95 : 1.0)
            .animation(.easeInOut(duration: 0.1), value: configuration.isPressed)
    }
}

// DNA-themed button style for special actions - using DNA helix green
struct GenomancerDNAButtonStyle: ButtonStyle {
    func makeBody(configuration: Configuration) -> some View {
        configuration.label
            .foregroundColor(.genomancerText)
            .padding()
            .background(
                RoundedRectangle(cornerRadius: 12)
                    .fill(Color.genomancerDNAHelix)
                    .opacity(configuration.isPressed ? 0.8 : 1.0)
            )
            .scaleEffect(configuration.isPressed ? 0.95 : 1.0)
            .animation(.easeInOut(duration: 0.1), value: configuration.isPressed)
    }
}

// Custom chip style for filter chips
struct GenomancerChipStyle: ViewModifier {
    let isSelected: Bool
    
    func body(content: Content) -> some View {
        content
            .foregroundColor(.genomancerText)
            .padding(.horizontal, 12)
            .padding(.vertical, 6)
            .background(
                RoundedRectangle(cornerRadius: 20)
                    .fill(isSelected ? Color.genomancerSecondary : Color.genomancerSecondaryBackground)
            )
    }
}

// Surface card style
struct GenomancerCardStyle: ViewModifier {
    func body(content: Content) -> some View {
        content
            .background(
                RoundedRectangle(cornerRadius: 16)
                    .fill(Color.genomancerCardBackground)
                    .shadow(color: .black.opacity(0.1), radius: 4, x: 0, y: 2)
            )
    }
}

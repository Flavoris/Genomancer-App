import SwiftUI
import DigestCore

struct GelView: View {
    let fragments: [Fragment]
    @State private var gelPercent: Double = 2.0
    
    var body: some View {
        VStack(spacing: 16) {
            // Gel percentage controls
            VStack(spacing: 12) {
                HStack {
                    Text("Gel %:")
                        .font(.system(size: 15, weight: .medium))
                        .foregroundColor(.white.opacity(0.8))
                    Spacer()
                    Text(String(format: "%.1f%%", gelPercent))
                        .font(.system(size: 15, weight: .semibold, design: .rounded))
                        .foregroundColor(.white)
                        .monospacedDigit()
                }
                .padding(.horizontal, 20)
                
                // Slider
                Slider(value: $gelPercent, in: 0.8...3.0, step: 0.1)
                    .tint(.blue)
                    .padding(.horizontal, 20)
                    .accessibilityLabel("Gel percentage")
                    .accessibilityValue(String(format: "%.1f percent", gelPercent))
                
                // Quick preset buttons
                HStack(spacing: 12) {
                    ForEach([1.0, 1.5, 2.0, 2.5, 3.0], id: \.self) { preset in
                        Button(action: {
                            gelPercent = preset
                        }) {
                            Text(String(format: "%.1f%%", preset))
                                .font(.system(size: 13, weight: .medium, design: .rounded))
                                .foregroundColor(gelPercent == preset ? .black : .white.opacity(0.7))
                                .frame(maxWidth: .infinity)
                                .padding(.vertical, 8)
                                .background(
                                    RoundedRectangle(cornerRadius: 8)
                                        .fill(gelPercent == preset ? Color.white : Color.white.opacity(0.1))
                                )
                        }
                    }
                }
                .padding(.horizontal, 20)
            }
            .padding(.vertical, 12)
            .background(Color.black.opacity(0.2))
            
            // Gel render view
            GelRenderView(fragments: fragments, gelPercent: gelPercent, style: .photo)
        }
        .background(Color.genomancerBackground)
        .navigationTitle("Gel Electrophoresis")
        .navigationBarTitleDisplayMode(.inline)
        .toolbarColorScheme(.dark, for: .navigationBar)
    }
}


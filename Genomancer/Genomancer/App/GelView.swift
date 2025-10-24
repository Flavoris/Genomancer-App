import SwiftUI
import DigestCore

struct GelView: View {
    let fragments: [Fragment]
    @State private var gelPercent: Double = 2.0
    
    var body: some View {
        HStack(spacing: 0) {
            // Gel render view
            GelRenderView(fragments: fragments, gelPercent: gelPercent, style: .photo)
            
            // Gel percentage controls on the right
            VStack(spacing: 12) {
                // Percentage display
                VStack(spacing: 4) {
                    Text("Gel %")
                        .font(.system(size: 13, weight: .medium))
                        .foregroundColor(.white.opacity(0.6))
                    Text(String(format: "%.1f%%", gelPercent))
                        .font(.system(size: 18, weight: .semibold, design: .rounded))
                        .foregroundColor(.white)
                        .monospacedDigit()
                }
                .padding(.top, 8)
                
                // Slider
                VStack {
                    Slider(value: $gelPercent, in: 0.5...3.0, step: 0.1)
                        .rotationEffect(.degrees(-90))
                        .frame(width: 200)
                        .tint(.blue)
                        .accessibilityLabel("Gel percentage")
                        .accessibilityValue(String(format: "%.1f percent", gelPercent))
                }
                .frame(width: 40, height: 200)
                .padding(.vertical, 16)
                
                // Quick preset buttons (vertical)
                VStack(spacing: 8) {
                    ForEach([3.0, 2.5, 2.0, 1.5, 1.0, 0.5], id: \.self) { preset in
                        Button(action: {
                            gelPercent = preset
                        }) {
                            Text(String(format: "%.1f%%", preset))
                                .font(.system(size: 12, weight: .medium, design: .rounded))
                                .foregroundColor(gelPercent == preset ? .black : .white.opacity(0.7))
                                .frame(width: 50, height: 32)
                                .background(
                                    RoundedRectangle(cornerRadius: 6)
                                        .fill(gelPercent == preset ? Color.white : Color.white.opacity(0.1))
                                )
                        }
                    }
                }
                .padding(.bottom, 8)
                
                Spacer()
            }
            .frame(width: 80)
        }
        .background(Color.genomancerBackground)
        .navigationTitle("Gel Electrophoresis")
        .navigationBarTitleDisplayMode(.inline)
        .toolbarColorScheme(.dark, for: .navigationBar)
    }
}


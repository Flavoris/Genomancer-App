import SwiftUI
import DigestCore

struct LigationCompatibilityView: View {
    let fragments: [Fragment]
    @State private var includeBlunt: Bool = false
    @State private var minOverhang: Int = 1
    @State private var requireDirectional: Bool = false
    @State private var selectedPair: CompatibilityResult?
    
    private var compatibilityResults: [CompatibilityResult] {
        calculateCompatibility(
            fragments: fragments,
            includeBlunt: includeBlunt,
            minOverhang: minOverhang,
            requireDirectional: requireDirectional
        )
    }
    
    var body: some View {
        VStack(spacing: 0) {
            // Options panel
            optionsPanel
            
            Divider()
                .background(Color.genomancerAccent.opacity(0.3))
            
            // Content area
            if compatibilityResults.isEmpty {
                emptyStateView
            } else {
                pairListView
            }
        }
        .background(Color.genomancerBackground)
        .navigationTitle("Ligation Compatibility")
        .toolbarColorScheme(.dark, for: .navigationBar)
        .sheet(item: $selectedPair) { pair in
            CompatibilityDetailSheet(result: pair, onDismiss: {
                selectedPair = nil
            })
        }
    }
    
    // MARK: - Options Panel
    
    private var optionsPanel: some View {
        VStack(spacing: 12) {
            // Options
            VStack(spacing: 8) {
                Toggle(isOn: $includeBlunt) {
                    Text("Include blunt-blunt pairs")
                        .font(.subheadline)
                        .foregroundColor(.genomancerText)
                }
                .tint(.genomancerAccent)
                
                Toggle(isOn: $requireDirectional) {
                    Text("Only directional pairs")
                        .font(.subheadline)
                        .foregroundColor(.genomancerText)
                }
                .tint(.genomancerAccent)
                
                HStack {
                    Text("Min overhang length:")
                        .font(.subheadline)
                        .foregroundColor(.genomancerText)
                    
                    Spacer()
                    
                    Stepper("\(minOverhang) bp", value: $minOverhang, in: 1...10)
                        .font(.subheadline)
                        .foregroundColor(.genomancerAccent)
                }
            }
            .padding(.horizontal)
            
            // Results summary
            HStack {
                Text("\(compatibilityResults.count) compatible pair(s)")
                    .font(.caption)
                    .foregroundColor(.genomancerSecondaryText)
                
                Spacer()
                
                if !compatibilityResults.isEmpty {
                    let directionalCount = compatibilityResults.filter { $0.directional }.count
                    Text("\(directionalCount) directional")
                        .font(.caption)
                        .foregroundColor(.genomancerAccent)
                }
            }
            .padding(.horizontal)
        }
        .padding(.vertical, 12)
        .background(Color.genomancerSecondaryBackground)
    }
    
    // MARK: - Empty State
    
    private var emptyStateView: some View {
        VStack(spacing: 16) {
            Image(systemName: "link.circle")
                .font(.system(size: 60))
                .foregroundColor(.genomancerSecondaryText)
            
            Text("No Compatible Pairs")
                .font(.title2)
                .fontWeight(.semibold)
                .foregroundColor(.genomancerText)
            
            Text("Try adjusting the options above to find compatible fragment ends.")
                .font(.subheadline)
                .foregroundColor(.genomancerSecondaryText)
                .multilineTextAlignment(.center)
                .padding(.horizontal, 40)
        }
        .frame(maxWidth: .infinity, maxHeight: .infinity)
    }
    
    // MARK: - Pair List View
    
    private var pairListView: some View {
        List(Array(compatibilityResults.enumerated()), id: \.offset) { index, result in
            Button(action: {
                selectedPair = result
            }) {
                CompatibilityPairRow(result: result, index: index)
            }
            .listRowInsets(EdgeInsets(top: 8, leading: 16, bottom: 8, trailing: 16))
            .listRowBackground(Color.clear)
            .listRowSeparator(.hidden)
        }
        .scrollContentBackground(.hidden)
        .background(Color.genomancerBackground)
        .listStyle(.plain)
    }
    
}

// MARK: - Compatibility Pair Row

struct CompatibilityPairRow: View {
    let result: CompatibilityResult
    let index: Int
    
    private var strengthColor: Color {
        if result.strengthScore >= 0.75 {
            return .green
        } else if result.strengthScore >= 0.5 {
            return .genomancerAccent
        } else {
            return .orange
        }
    }
    
    var body: some View {
        VStack(alignment: .leading, spacing: 16) {
            // Header with strength indicator
            HStack {
                Text("Pair #\(index + 1)")
                    .font(.headline)
                    .foregroundColor(.genomancerAccent)
                
                Spacer()
                
                // Strength badge
                HStack(spacing: 4) {
                    Circle()
                        .fill(strengthColor)
                        .frame(width: 8, height: 8)
                    Text("\(Int(result.strengthScore * 100))%")
                        .font(.caption)
                        .fontWeight(.semibold)
                        .foregroundColor(.white)
                }
                .padding(.horizontal, 10)
                .padding(.vertical, 6)
                .background(Color.genomancerBackground.opacity(0.6))
                .cornerRadius(8)
                
                if result.directional {
                    Label("Directional", systemImage: "arrow.right")
                        .font(.caption)
                        .fontWeight(.semibold)
                        .foregroundColor(.green)
                        .padding(.horizontal, 8)
                        .padding(.vertical, 4)
                        .background(Color.genomancerBackground.opacity(0.6))
                        .cornerRadius(6)
                } else {
                    Label("Palindromic", systemImage: "arrow.left.arrow.right")
                        .font(.caption)
                        .fontWeight(.semibold)
                        .foregroundColor(.orange)
                        .padding(.horizontal, 8)
                        .padding(.vertical, 4)
                        .background(Color.genomancerBackground.opacity(0.6))
                        .cornerRadius(6)
                }
            }
            .padding(.bottom, 4)
            
            // End A
            endInfoView(
                label: "End A",
                end: result.endA,
                gcPercent: result.gcPercentA,
                tm: result.tmA
            )
            
            // Compatibility indicator with divider
            VStack(spacing: 0) {
                Divider()
                    .background(Color.genomancerAccent.opacity(0.3))
                
                HStack {
                    Spacer()
                    Image(systemName: "arrow.up.arrow.down.circle.fill")
                        .font(.title3)
                        .foregroundColor(.genomancerAccent)
                    Spacer()
                }
                .padding(.vertical, 8)
                
                Divider()
                    .background(Color.genomancerAccent.opacity(0.3))
            }
            .padding(.vertical, 8)
            
            // End B
            endInfoView(
                label: "End B",
                end: result.endB,
                gcPercent: result.gcPercentB,
                tm: result.tmB
            )
            
            // Note with better visual separation
            VStack(alignment: .leading, spacing: 4) {
                HStack {
                    Image(systemName: "info.circle.fill")
                        .font(.caption)
                        .foregroundColor(.genomancerAccent.opacity(0.7))
                    Text("Details")
                        .font(.caption2)
                        .fontWeight(.semibold)
                        .foregroundColor(.genomancerSecondaryText)
                        .textCase(.uppercase)
                }
                
                Text(result.note)
                    .font(.caption)
                    .foregroundColor(.genomancerText)
            }
            .frame(maxWidth: .infinity, alignment: .leading)
            .padding(.horizontal, 12)
            .padding(.vertical, 10)
            .background(Color.genomancerBackground.opacity(0.8))
            .cornerRadius(8)
            .overlay(
                RoundedRectangle(cornerRadius: 8)
                    .stroke(Color.genomancerAccent.opacity(0.2), lineWidth: 1)
            )
            .padding(.top, 4)
        }
        .padding(16)
        .background(Color(UIColor.secondarySystemBackground))
        .cornerRadius(12)
        .shadow(color: Color.black.opacity(0.2), radius: 4, x: 0, y: 2)
    }
    
    private func endInfoView(label: String, end: FragmentEndInfo,
                            gcPercent: Double, tm: Double) -> some View {
        VStack(alignment: .leading, spacing: 8) {
            // Label badge
            HStack {
                Text(label)
                    .font(.caption)
                    .fontWeight(.semibold)
                    .foregroundColor(.white)
                    .textCase(.uppercase)
                    .padding(.horizontal, 8)
                    .padding(.vertical, 3)
                    .background(Color.genomancerAccent.opacity(0.7))
                    .cornerRadius(4)
                Spacer()
            }
            
            // Fragment info
            HStack {
                Text("Fragment \(end.fragmentId + 1)")
                    .font(.subheadline)
                    .fontWeight(.semibold)
                    .foregroundColor(.genomancerText)
                
                Text("(\(end.polarity) end)")
                    .font(.caption)
                    .foregroundColor(.genomancerSecondaryText)
                
                if let enzyme = end.endInfo.sourceEnzyme {
                    Text("•")
                        .foregroundColor(.genomancerSecondaryText)
                    Text(enzyme)
                        .font(.caption)
                        .fontWeight(.semibold)
                        .foregroundColor(.white)
                        .padding(.horizontal, 6)
                        .padding(.vertical, 2)
                        .background(Color.genomancerAccent.opacity(0.7))
                        .cornerRadius(4)
                }
            }
            
            // Overhang sequence and properties
            if let seq = end.endInfo.overhangSeq5to3, !seq.isEmpty {
                VStack(alignment: .leading, spacing: 6) {
                    Text("5'-\(seq)-3'")
                        .font(.system(.subheadline, design: .monospaced))
                        .fontWeight(.medium)
                        .foregroundColor(.genomancerText)
                        .padding(.horizontal, 10)
                        .padding(.vertical, 6)
                        .background(Color.genomancerBackground.opacity(0.6))
                        .cornerRadius(6)
                    
                    HStack(spacing: 12) {
                        Label {
                            Text("\(gcPercent, specifier: "%.1f")%")
                                .font(.caption)
                                .foregroundColor(.genomancerText)
                        } icon: {
                            Text("GC:")
                                .font(.caption2)
                                .foregroundColor(.genomancerSecondaryText)
                        }
                        
                        Label {
                            Text("\(tm, specifier: "%.0f")°C")
                                .font(.caption)
                                .foregroundColor(.genomancerText)
                        } icon: {
                            Text("Tm≈")
                                .font(.caption2)
                                .foregroundColor(.genomancerSecondaryText)
                        }
                    }
                }
            }
        }
        .padding(12)
        .frame(maxWidth: .infinity, alignment: .leading)
        .background(Color.genomancerBackground.opacity(0.4))
        .cornerRadius(10)
        .overlay(
            RoundedRectangle(cornerRadius: 10)
                .stroke(Color.genomancerAccent.opacity(0.15), lineWidth: 1)
        )
    }
}

// MARK: - Compatibility Detail Sheet

struct CompatibilityDetailSheet: View {
    let result: CompatibilityResult
    let onDismiss: () -> Void
    
    private var strengthColor: Color {
        if result.strengthScore >= 0.75 {
            return .green
        } else if result.strengthScore >= 0.5 {
            return .genomancerAccent
        } else {
            return .orange
        }
    }
    
    var body: some View {
        NavigationView {
            ScrollView {
                VStack(spacing: 24) {
                    // Strength indicator
                    VStack(spacing: 12) {
                        Text("Compatibility Strength")
                            .font(.headline)
                            .foregroundColor(.genomancerText)
                        
                        ZStack(alignment: .leading) {
                            RoundedRectangle(cornerRadius: 8)
                                .fill(Color.genomancerSecondaryBackground)
                                .frame(height: 40)
                            
                            RoundedRectangle(cornerRadius: 8)
                                .fill(strengthColor)
                                .frame(width: CGFloat(result.strengthScore) * (UIScreen.main.bounds.width - 80), height: 40)
                            
                            HStack {
                                Text("\(Int(result.strengthScore * 100))%")
                                    .font(.headline)
                                    .fontWeight(.bold)
                                    .foregroundColor(.white)
                                    .padding(.leading, 12)
                                
                                Spacer()
                                
                                Text(result.strengthLabel)
                                    .font(.subheadline)
                                    .fontWeight(.semibold)
                                    .foregroundColor(result.strengthScore > 0.3 ? .white : .genomancerText)
                                    .padding(.trailing, 12)
                            }
                        }
                    }
                    .padding(.horizontal)
                    
                    // Compatibility type
                    VStack(spacing: 8) {
                        HStack {
                            if result.directional {
                                Label("Directional Cloning", systemImage: "arrow.right.circle.fill")
                                    .foregroundColor(.green)
                            } else {
                                Label("Non-Directional (Palindromic)", systemImage: "arrow.left.arrow.right.circle.fill")
                                    .foregroundColor(.orange)
                            }
                            Spacer()
                        }
                        .font(.subheadline)
                        .padding()
                        .background(Color.genomancerSecondaryBackground)
                        .cornerRadius(12)
                        
                        Text(result.note)
                            .font(.caption)
                            .foregroundColor(.genomancerSecondaryText)
                            .frame(maxWidth: .infinity, alignment: .leading)
                    }
                    .padding(.horizontal)
                    
                    Divider()
                        .background(Color.genomancerAccent.opacity(0.3))
                    
                    // End A details
                    endDetailCard(
                        title: "End A",
                        end: result.endA,
                        gcPercent: result.gcPercentA,
                        tm: result.tmA
                    )
                    
                    // Compatibility symbol
                    HStack {
                        Spacer()
                        VStack(spacing: 4) {
                            Image(systemName: "link.circle.fill")
                                .font(.system(size: 32))
                                .foregroundColor(strengthColor)
                            Text("Compatible")
                                .font(.caption)
                                .foregroundColor(.genomancerSecondaryText)
                        }
                        Spacer()
                    }
                    
                    // End B details
                    endDetailCard(
                        title: "End B",
                        end: result.endB,
                        gcPercent: result.gcPercentB,
                        tm: result.tmB
                    )
                    
                    Divider()
                        .background(Color.genomancerAccent.opacity(0.3))
                    
                    // Analysis breakdown
                    VStack(alignment: .leading, spacing: 16) {
                        Text("Strength Analysis")
                            .font(.headline)
                            .foregroundColor(.genomancerText)
                        
                        analysisRow(
                            label: "GC Content",
                            value: String(format: "%.1f%%", (result.gcPercentA + result.gcPercentB) / 2.0),
                            optimal: "~50%",
                            score: 1.0 - abs((result.gcPercentA + result.gcPercentB) / 2.0 - 50.0) / 50.0
                        )
                        
                        analysisRow(
                            label: "Melting Temp",
                            value: String(format: "~%.0f°C", (result.tmA + result.tmB) / 2.0),
                            optimal: ">30°C",
                            score: min((result.tmA + result.tmB) / 2.0 / 50.0, 1.0)
                        )
                        
                        analysisRow(
                            label: "Overhang Length",
                            value: "\(result.endA.endInfo.overhangLen) bp",
                            optimal: "4-6 bp",
                            score: min(Double(result.endA.endInfo.overhangLen) / 6.0, 1.0)
                        )
                        
                        analysisRow(
                            label: "Directionality",
                            value: result.directional ? "Yes" : "No",
                            optimal: "Yes",
                            score: result.directional ? 1.0 : 0.0
                        )
                    }
                    .padding()
                    .background(Color.genomancerSecondaryBackground)
                    .cornerRadius(12)
                    .padding(.horizontal)
                }
                .padding(.vertical)
            }
            .background(Color.genomancerBackground)
            .navigationTitle("Compatibility Details")
            .navigationBarTitleDisplayMode(.inline)
            .toolbarColorScheme(.dark, for: .navigationBar)
            .toolbar {
                ToolbarItem(placement: .navigationBarTrailing) {
                    Button("Done") {
                        onDismiss()
                    }
                    .foregroundColor(.genomancerAccent)
                }
            }
        }
    }
    
    private func endDetailCard(title: String, end: FragmentEndInfo,
                               gcPercent: Double, tm: Double) -> some View {
        VStack(alignment: .leading, spacing: 12) {
            Text(title)
                .font(.headline)
                .foregroundColor(.genomancerAccent)
            
            VStack(alignment: .leading, spacing: 8) {
                HStack {
                    Text("Fragment \(end.fragmentId + 1)")
                        .font(.title3)
                        .fontWeight(.semibold)
                        .foregroundColor(.genomancerText)
                    
                    Text("(\(end.polarity) end)")
                        .font(.subheadline)
                        .foregroundColor(.genomancerSecondaryText)
                }
                
                if let enzyme = end.endInfo.sourceEnzyme {
                    HStack {
                        Text("Enzyme:")
                            .font(.subheadline)
                            .foregroundColor(.genomancerSecondaryText)
                        Text(enzyme)
                            .font(.subheadline)
                            .fontWeight(.semibold)
                            .foregroundColor(.genomancerAccent)
                    }
                }
                
                HStack {
                    Text("Position:")
                        .font(.subheadline)
                        .foregroundColor(.genomancerSecondaryText)
                    Text("\(end.position) bp")
                        .font(.subheadline)
                        .foregroundColor(.genomancerText)
                }
                
                VStack(alignment: .leading, spacing: 4) {
                    Text("Overhang:")
                        .font(.subheadline)
                        .foregroundColor(.genomancerSecondaryText)
                    
                    if let seq = end.endInfo.overhangSeq5to3, !seq.isEmpty {
                        Text("5'-\(seq)-3'")
                            .font(.system(.title3, design: .monospaced))
                            .fontWeight(.bold)
                            .foregroundColor(.genomancerText)
                            .padding(.vertical, 8)
                            .padding(.horizontal, 12)
                            .background(Color.genomancerBackground)
                            .cornerRadius(8)
                        
                        HStack(spacing: 20) {
                            VStack(alignment: .leading, spacing: 2) {
                                Text("GC Content")
                                    .font(.caption)
                                    .foregroundColor(.genomancerSecondaryText)
                                Text("\(gcPercent, specifier: "%.1f")%")
                                    .font(.subheadline)
                                    .fontWeight(.semibold)
                                    .foregroundColor(.genomancerText)
                            }
                            
                            VStack(alignment: .leading, spacing: 2) {
                                Text("Melting Temp")
                                    .font(.caption)
                                    .foregroundColor(.genomancerSecondaryText)
                                Text("~\(tm, specifier: "%.0f")°C")
                                    .font(.subheadline)
                                    .fontWeight(.semibold)
                                    .foregroundColor(.genomancerText)
                            }
                            
                            VStack(alignment: .leading, spacing: 2) {
                                Text("Type")
                                    .font(.caption)
                                    .foregroundColor(.genomancerSecondaryText)
                                Text(end.endInfo.overhangType == .fivePrime ? "5' overhang" : "3' overhang")
                                    .font(.subheadline)
                                    .fontWeight(.semibold)
                                    .foregroundColor(.genomancerText)
                            }
                        }
                    } else {
                        Text("Blunt")
                            .font(.subheadline)
                            .foregroundColor(.genomancerText)
                    }
                }
            }
        }
        .padding()
        .background(Color.genomancerSecondaryBackground)
        .cornerRadius(12)
        .padding(.horizontal)
    }
    
    private func analysisRow(label: String, value: String, optimal: String, score: Double) -> some View {
        VStack(spacing: 8) {
            HStack {
                Text(label)
                    .font(.subheadline)
                    .foregroundColor(.genomancerText)
                Spacer()
                Text(value)
                    .font(.subheadline)
                    .fontWeight(.semibold)
                    .foregroundColor(.genomancerAccent)
            }
            
            HStack {
                Text("Optimal: \(optimal)")
                    .font(.caption)
                    .foregroundColor(.genomancerSecondaryText)
                
                Spacer()
                
                // Score bar
                ZStack(alignment: .leading) {
                    RoundedRectangle(cornerRadius: 4)
                        .fill(Color.genomancerBackground)
                        .frame(width: 100, height: 8)
                    
                    RoundedRectangle(cornerRadius: 4)
                        .fill(score >= 0.75 ? Color.green : (score >= 0.5 ? Color.genomancerAccent : Color.orange))
                        .frame(width: CGFloat(score) * 100, height: 8)
                }
            }
        }
    }
}

// MARK: - Preview

struct LigationCompatibilityView_Previews: PreviewProvider {
    static var previews: some View {
        NavigationView {
            LigationCompatibilityView(fragments: [
                Fragment(
                    start: 0,
                    end: 100,
                    length: 100,
                    leftEnd: EndInfo(overhangType: .fivePrime, overhangSeq5to3: "AATT", sourceEnzyme: "EcoRI", overhangLen: 4),
                    rightEnd: EndInfo(overhangType: .fivePrime, overhangSeq5to3: "GATC", sourceEnzyme: "BamHI", overhangLen: 4),
                    sequence: "ATCGATCG"
                ),
                Fragment(
                    start: 100,
                    end: 200,
                    length: 100,
                    leftEnd: EndInfo(overhangType: .fivePrime, overhangSeq5to3: "GATC", sourceEnzyme: "BamHI", overhangLen: 4),
                    rightEnd: EndInfo(overhangType: .fivePrime, overhangSeq5to3: "AATT", sourceEnzyme: "EcoRI", overhangLen: 4),
                    sequence: "GCTAGCTA"
                )
            ])
        }
    }
}


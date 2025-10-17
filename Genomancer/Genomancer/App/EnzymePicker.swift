import SwiftUI
import DigestCore

struct EnzymePicker: View {
    let all: [Enzyme]
    @Binding var selected: Set<Enzyme>
    @State private var query = ""
    @State private var overhangFilter: OverhangType? = nil

    var body: some View {
        VStack(spacing: 0) {
            // Filter chips
            ScrollView(.horizontal, showsIndicators: false) {
                HStack(spacing: 8) {
                    FilterChip(
                        label: "All",
                        isSelected: overhangFilter == nil,
                        action: { overhangFilter = nil }
                    )
                    FilterChip(
                        label: "Blunt",
                        isSelected: overhangFilter == .blunt,
                        action: { overhangFilter = .blunt }
                    )
                    FilterChip(
                        label: "5' Overhang",
                        isSelected: overhangFilter == .fivePrime,
                        action: { overhangFilter = .fivePrime }
                    )
                    FilterChip(
                        label: "3' Overhang",
                        isSelected: overhangFilter == .threePrime,
                        action: { overhangFilter = .threePrime }
                    )
                }
                .padding(.horizontal)
                .padding(.vertical, 8)
            }
            .background(Color.genomancerCardBackground)
            
            List(filtered, id: \.self) { e in
            HStack {
                VStack(alignment: .leading) {
                    Text(e.name)
                        .font(.headline)
                        .dynamicTypeSize(...DynamicTypeSize.accessibility3)
                    Text(e.site)
                        .font(.caption)
                        .foregroundColor(.genomancerSecondaryText)
                        .dynamicTypeSize(...DynamicTypeSize.accessibility3)
                }
                Spacer()
                if selected.contains(e) { 
                    Image(systemName: "checkmark")
                        .accessibilityHidden(true)
                }
            }
            .contentShape(Rectangle())
            .onTapGesture {
                if selected.contains(e) { selected.remove(e) } else { selected.insert(e) }
            }
            .accessibilityElement(children: .combine)
            .accessibilityLabel(e.name)
            .accessibilityValue("Recognition site: \(e.site), \(overhangLabel(e.overhangType))")
            .accessibilityHint(selected.contains(e) ? "Selected. Double tap to deselect" : "Not selected. Double tap to select")
            .accessibilityAddTraits(selected.contains(e) ? [.isSelected] : [])
        }
        .searchable(text: $query)
        }
        .background(Color.genomancerBackground)
        .navigationTitle("Enzymes")
    }

    var filtered: [Enzyme] {
        var enzymes = all
        
        // Apply overhang filter
        if let filter = overhangFilter {
            enzymes = enzymes.filter { $0.overhangType == filter }
        }
        
        // Apply search query
        if !query.isEmpty {
            let q = query.lowercased()
            enzymes = enzymes.filter { $0.name.lowercased().contains(q) || $0.site.lowercased().contains(q) }
        }
        
        return enzymes
    }
    
    private func overhangLabel(_ type: OverhangType?) -> String {
        guard let type = type else { return "unknown overhang type" }
        switch type {
        case .blunt:
            return "blunt ends"
        case .fivePrime:
            return "5' overhang"
        case .threePrime:
            return "3' overhang"
        case .unknown:
            return "unknown overhang type"
        }
    }
}

struct FilterChip: View {
    let label: String
    let isSelected: Bool
    let action: () -> Void
    
    var body: some View {
        Button(action: action) {
            Text(label)
                .font(.subheadline)
                .padding(.horizontal, 12)
                .padding(.vertical, 6)
                .background(isSelected ? Color.genomancerSecondary : Color.genomancerSecondaryBackground)
                .foregroundColor(.genomancerText)
                .cornerRadius(16)
        }
    }
}


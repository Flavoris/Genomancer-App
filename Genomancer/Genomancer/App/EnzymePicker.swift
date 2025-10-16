import SwiftUI
import DigestCore

struct EnzymePicker: View {
    let all: [Enzyme]
    @Binding var selected: Set<Enzyme>
    @State private var query = ""

    var body: some View {
        List(filtered, id: \.self) { e in
            HStack {
                VStack(alignment: .leading) {
                    Text(e.name).font(.headline)
                    Text(e.site).font(.caption).foregroundStyle(.secondary)
                }
                Spacer()
                if selected.contains(e) { Image(systemName: "checkmark") }
            }
            .contentShape(Rectangle())
            .onTapGesture {
                if selected.contains(e) { selected.remove(e) } else { selected.insert(e) }
            }
        }
        .searchable(text: $query)
        .navigationTitle("Enzymes")
    }

    var filtered: [Enzyme] {
        guard !query.isEmpty else { return all }
        let q = query.lowercased()
        return all.filter { $0.name.lowercased().contains(q) || $0.site.lowercased().contains(q) }
    }
}


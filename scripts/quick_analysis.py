# Quick analysis of sequence lengths
with open('data/extracted_sequences.fasta', 'r') as f:
    lines = f.readlines()

lengths = []
for line in lines:
    if line.startswith('>'):
        # Extract coordinates (e.g., chr12:53676107-53676353)
        if ':' in line and '-' in line:
            coords = line.strip().split(':')[1]
            if '-' in coords:
                start, end = coords.split('-')
                length = int(end) - int(start) + 1
                lengths.append(length)

print(f'Analysis of {len(lengths)} sequences:')
print(f'Min length: {min(lengths)} bp')
print(f'Max length: {max(lengths)} bp')
print(f'Average length: {sum(lengths)/len(lengths):.1f} bp')

# Length distribution
short = sum(1 for l in lengths if l <= 100)
medium = sum(1 for l in lengths if 101 <= l <= 200)
long = sum(1 for l in lengths if 201 <= l <= 300)
very_long = sum(1 for l in lengths if l > 300)

print(f'Length distribution:')
print(f'  <= 100 bp: {short:,} sequences ({short/len(lengths)*100:.1f}%)')
print(f'  101-200 bp: {medium:,} sequences ({medium/len(lengths)*100:.1f}%)')
print(f'  201-300 bp: {long:,} sequences ({long/len(lengths)*100:.1f}%)')
print(f'  > 300 bp: {very_long:,} sequences ({very_long/len(lengths)*100:.1f}%)')

print(f'\nWith current max_length=100: Only {short:,} sequences would pass ({short/len(lengths)*100:.1f}%)')
print(f'With max_length=200: {short+medium:,} sequences would pass ({(short+medium)/len(lengths)*100:.1f}%)')
print(f'With max_length=300: {short+medium+long:,} sequences would pass ({(short+medium+long)/len(lengths)*100:.1f}%)')

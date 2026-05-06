import re

# Read the file
with open('paper.tex', 'rb') as f:
    content = f.read()

# The pattern is a literal backslash-dollar-whitespace-ext-AA
# We'll search for any variant of: \$ followed by whitespace and ext{AA}$$

# Convert to string for processing
content_str = content.decode('utf-8', errors='replace')

# Replace all variations of the broken pattern
content_str = content_str.replace(r'\$ ext{AA}$$', r'$\text{\AA}$')
content_str = content_str.replace(r'\$	ext{AA}$$', r'$\text{\AA}$')
content_str = content_str.replace(r'\$ ext{AA}$', r'$\text{\AA}$')
content_str = content_str.replace(r'\$	ext{AA}$', r'$\text{\AA}$')

# Also fix the cases with explicit whitespace shown
content_str = re.sub(r'\\$\s+ext\{AA\}\$\$', r'$\text{\AA}$', content_str)
content_str = re.sub(r'\\$\s+ext\{AA\}\$', r'$\text{\AA}$', content_str)

with open('paper.tex', 'w', encoding='utf-8') as f:
    f.write(content_str)

print("Fixed!")

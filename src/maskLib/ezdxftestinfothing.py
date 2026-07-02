import ezdxf
from ezdxf.addons import text2path
from ezdxf.tools.fonts import FontFace
from dxfwrite import DXFEngine as dxf

# Define the text and font
text = "p"
font_face = FontFace(family='Arial')

# Convert text to paths
paths = text2path.make_paths_from_str(text, font=font_face, size=10)

# Get the starting point of the first path (which should be the "p")
start_point = paths[0].start

print(f"Starting point of 'p': {start_point}")

# Create a new DXF document
doc = ezdxf.new('R2010')
msp = doc.modelspace()

# Create a block for the text
block_name = "TEXT_BLOCK"
text_block = dxf.block(block_name)

# Add the text paths to the block using the same polyline command as in maskLib.py
for path in paths:
    points = list(path.flattening(0.01))
    # Ensure the path is closed by adding the starting point at the end
    if points[0] != points[-1]:
        points.append(points[0])
    text_block.add(dxf.polyline(points, layer='TEXT', flags=1))  # flags=1 ensures the polyline is closed

# Add the block to the drawing
doc.blocks.add(text_block)

# Add the block reference to the modelspace
msp.add_blockref(block_name, insert=(0, 0))

# Add a circle at the starting point for easy identification
msp.add_circle(center=(start_point.x, start_point.y), radius=0.1, dxfattribs={'color': 1})

# Save the DXF document
doc.saveas('text_with_start_point.dxf')
print("DXF file created: text_with_start_point.dxf")
'''
project: sequencing
april 29 2025 - XQ - Leiden UMC
generate colour scheme (from R to python)
''' 

import colorsys
import numpy as np


def generate_scheme(base_color, n=5, light_range=(0.4, 0.1), dark_range=(0.1, 0.4)):
    """
    Generate Monochromatic Color Schemes

    Return one or more lists of lighter and darker shades based on one or more base colors.
    If multiple base colors are supplied, a dictionary of color schemes will be returned.

    Parameters:
    -----------
    base_color : str or list
        A hex color string or a list of hex strings (e.g., "#7FC97F" or ["#7FC97F", "#BEAED4"]).
    n : int, optional
        Number of colors in each scheme (default is 5).
    light_range : tuple, optional
        A tuple of length 2 specifying the range of lightening (default = (0.4, 0.1)).
    dark_range : tuple, optional
        A tuple of length 2 specifying the range of darkening (default = (0.1, 0.4)).

    Returns:
    --------
    list or dict
        A list of hex colors if one base color is provided, or a dictionary of such lists if multiple.

    Examples:
    --------
    >>> generate_scheme("#7FC97F")
    >>> generate_scheme(["#7FC97F", "#BEAED4"], n=7, light_range=(0.6, 0.2), dark_range=(0.2, 0.6))
    """
    def hex_to_rgb(hex_str):
        hex_str = hex_str.lstrip('#')
        return tuple(int(hex_str[i:i+2], 16) / 255.0 for i in (0, 2, 4))
    
    def rgb_to_hex(rgb):
        return '#{:02x}{:02x}{:02x}'.format(
            int(min(max(rgb[0] * 255, 0), 255)),
            int(min(max(rgb[1] * 255, 0), 255)),
            int(min(max(rgb[2] * 255, 0), 255))
        )
    
    def adjust_lightness(color_rgb, amount, lighten=True):
        h, l, s = colorsys.rgb_to_hls(*color_rgb)
        if lighten:
            l = min(l + amount, 1.0)
        else:
            l = max(l - amount, 0.0)
        r, g, b = colorsys.hls_to_rgb(h, l, s)
        return (r, g, b)
    
    def generate_single(clr):
        n_light = int(np.ceil(n / 2))
        n_dark = int(np.floor(n / 2))
        
        lighten_vals = np.linspace(light_range[0], light_range[1], n_light)
        darken_vals = np.linspace(dark_range[0], dark_range[1], n_dark)
        
        rgb = hex_to_rgb(clr)
        light_colors = [rgb_to_hex(adjust_lightness(rgb, a, lighten=True)) for a in lighten_vals]
        dark_colors = [rgb_to_hex(adjust_lightness(rgb, a, lighten=False)) for a in darken_vals]
        
        return light_colors + dark_colors
    
    if isinstance(base_color, str):
        return generate_single(base_color)
    else:
        schemes = {color: generate_single(color) for color in base_color}
        return schemes
# FlowForge

**FlowForge** is an independent research project focused on generating a large dataset of 250+ realistic winglet geometries across eight major categories. Including traditional, blended, raked, and spiroid designs. For use in CFD simulations and machine learning–based aerodynamic optimization.

## Usage
1. Run the generation script to create winglet geometries:
   ```bash
   python generate_winglets.py --output ./winglets
   
2. All models are automatically exported as .stl files.

3. Each category folder contains realistic geometric variations ready for CFD data generation

## Features

- Procedural generation of 250+ winglet models
- Includes 8 winglet categories with smooth, realistic variations
- Exports .stl files intended for CFD meshing
- Designed for ML and surrogate-based aerodynamic studies
- Fully self-contained, no OpenVSP or external modeling software required
- 
## Contributing

If you’d like to contribute improvements or additional geometry types:

1. Fork the repository
2. Create a new branch
3. Submit a pull request

Things that may still need improvement

- Smoother blending of winglet into airfoil
- better forms of the more complicated winglet types, like wing grids

## Licence

MIT License

Copyright (c) 2025 Advaith Reddy Voodem

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including, without limitation, the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT, OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OF OTHER DEALINGS IN THE
SOFTWARE.

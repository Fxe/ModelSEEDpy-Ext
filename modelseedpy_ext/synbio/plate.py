import logging

logger = logging.getLogger(__name__)


class PlateWell:

    def __init__(self, medium='?', culture='?', measurements=None):
        self.medium = medium
        self.culture = culture
        self.measurements = measurements
        if self.measurements is None:
            self.measurements = []

    def add_measurement(self, time, kind, value):
        """Add a measurement to this well"""
        self.measurements.append({
            'time': time,
            'kind': kind,
            'value': value
        })


class WellPlate96:

    def __init__(self, plate_name):
        """Initialize a 96-well plate with wells A1-H12"""

        self.wells = {}
        self.rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        self.cols = list(range(1, 13))  # 1-12
        self.plate_name = plate_name
        self.wells_style = {}  # display style
        
        # Create all 96 wells
        for row in self.rows:
            for col in self.cols:
                well_id = f"{row}{col}"
                self.wells[well_id] = PlateWell()
                self.wells_style[well_id] = {}

    def get_well(self, row, col):
        """
        Get a well by row and column.
        
        Args:
            row: Either string ('A'-'H') or int (0-7)
            col: Either int (1-12) or int (0-11) 
            
        Returns:
            PlateWell: The well at the specified position
        """
        # Convert row to letter if it's an integer
        if isinstance(row, int):
            if 0 <= row <= 7:
                row = self.rows[row]
            else:
                raise ValueError(f"Row index {row} out of range (0-7)")
        
        # Convert column to 1-12 range if it's 0-11
        if isinstance(col, int):
            if 0 <= col <= 11:
                col = col + 1
            elif 1 <= col <= 12:
                pass  # Already in correct range
            else:
                raise ValueError(f"Column index {col} out of range (0-11 or 1-12)")
        
        well_id = f"{row}{col}"
        if well_id not in self.wells:
            raise ValueError(f"Invalid well position: {well_id}")
            
        return self.wells[well_id]
    
    def set_well(self, row, col, medium='?', culture='?', measurements=[]):
        """
        Set or update a well with new parameters.
        
        Args:
            row: Row identifier ('A'-'H' or 0-7)
            col: Column identifier (1-12 or 0-11)
            medium: Medium type
            culture: Culture/microorganism
            measurements: List of measurements
        """
        well = self.get_well(row, col)
        well.medium = medium
        well.culture = culture
        well.measurements = measurements if measurements else []
        
    def get_well_id(self, row, col):
        """Get the well ID string (e.g., 'A1') for given row/col"""
        if isinstance(row, int):
            if 0 <= row <= 7:
                row = self.rows[row]
            else:
                raise ValueError(f"Row index {row} out of range (0-7)")
                
        if isinstance(col, int):
            if 0 <= col <= 11:
                col = col + 1
            elif 1 <= col <= 12:
                pass
            else:
                raise ValueError(f"Column index {col} out of range (0-11 or 1-12)")
                
        return f"{row}{col}"
    
    def __iter__(self):
        """Make the plate iterable - yields wells in row-major order (A1-A12, B1-B12, etc.)"""
        for row in self.rows:
            for col in self.cols:
                well_id = f"{row}{col}"
                yield self.wells[well_id]
                
    def get_row(self, row):
        """Get all wells in a specific row"""
        if isinstance(row, int):
            if 0 <= row <= 7:
                row = self.rows[row]
            else:
                raise ValueError(f"Row index {row} out of range (0-7)")
                
        return [self.wells[f"{row}{col}"] for col in self.cols]
    
    def get_column(self, col):
        """Get all wells in a specific column"""
        if isinstance(col, int):
            if 0 <= col <= 11:
                col = col + 1
            elif 1 <= col <= 12:
                pass
            else:
                raise ValueError(f"Column index {col} out of range (0-11 or 1-12)")
                
        return [self.wells[f"{row}{col}"] for row in self.rows]
    
    def _generate_colors(self, unique_values):
        """Generate distinct colors for different values"""
        colors = [
            '#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7',
            '#DDA0DD', '#98D8C8', '#F7DC6F', '#BB8FCE', '#85C1E9',
            '#F8C471', '#82E0AA', '#F1948A', '#85C1E9', '#D7DBDD',
            '#AED6F1', '#A9DFBF', '#F9E79F', '#D2B4DE', '#AEB6BF'
        ]
        
        color_map = {}
        for i, value in enumerate(unique_values):
            if value == '?':  # Default/empty wells
                color_map[value] = '#F8F9FA'  # Light gray
            else:
                color_map[value] = colors[i % len(colors)]
        
        return color_map

    def get_culture_list(self):
        """Get a list of all unique cultures in the plate"""
        return list(set([well.culture for well in self]))
    
    def get_medium_list(self):
        """Get a list of all unique media types in the plate"""
        return list(set([well.medium for well in self]))

    def _repr_html_(self, color_by='medium', box_size=50):
        """
        Generate HTML representation of the plate for Jupyter notebooks
        
        Args:
            color_by (str): 'medium' or 'culture' - what to color the wells by
            box_size (int): Size of each well box in pixels (default: 50)
        """
        # Collect unique values for coloring
        if color_by == 'medium':
            values = [str(x) for x in self.get_medium_list()]
            title_suffix = "Medium Types"
        elif color_by == 'culture':
            values = [str(x) for x in self.get_culture_list()]
            title_suffix = "Culture Types"
        else:
            raise ValueError("color_by must be 'medium' or 'culture'")
            
        unique_values = sorted(list(set(values)))
        color_map = self._generate_colors(unique_values)
        
        # Start building HTML
        html = f"""
        <div style="font-family: Arial, sans-serif; margin: 20px;">
            <h3 style="text-align: center; color: #333;">
                96-Well Plate: {self.plate_name} - Colored by {title_suffix}
            </h3>
            
            <!-- Legend -->
            <div style="margin-bottom: 20px; text-align: center;">
                <strong>Legend:</strong><br>
                <div style="display: inline-block; text-align: left;">
        """

        # Add legend items
        for value in unique_values:
            color = color_map[value]
            display_name = value if value != '?' else 'Empty/Unknown'
            html += f"""
                    <span style="display: inline-block; margin: 2px 10px;">
                        <span style="width: 15px; height: 15px; background-color: {color}; 
                                   display: inline-block; border: 1px solid #ccc; margin-right: 5px;"></span>
                        {display_name}
                    </span>
            """
        
        html += """
                </div>
            </div>
            
            <!-- Plate Layout -->
            <table style="border-collapse: collapse; margin: 0 auto; border: 2px solid #333;">
                <thead>
                    <tr>
                        <th style="background-color: #f0f0f0; border: 1px solid #ccc; padding: 8px; font-weight: bold;"></th>
        """
        
        # Column headers (1-12)
        for col in self.cols:
            html += f'<th style="background-color: #f0f0f0; border: 1px solid #ccc; padding: 8px; font-weight: bold; text-align: center;">{col}</th>'
        
        html += """
                    </tr>
                </thead>
                <tbody>
        """
        
        # Generate plate rows
        for row in self.rows:
            html += f"""
                    <tr>
                        <td style="background-color: #f0f0f0; border: 1px solid #ccc; padding: 8px; font-weight: bold; text-align: center;">{row}</td>
            """
            
            for col in self.cols:
                well_id = f"{row}{col}"
                well = self.wells[well_id]
                well_style = self.wells_style[well_id]
                if color_by == 'medium':
                    value = well.medium
                else:
                    value = well.culture
                value = str(value)
                color = color_map[value]
                display_value = value if value != '?' else ''
                # Truncate long values for display
                if len(display_value) > 8:
                    display_value = display_value[:6] + '...'
                border_color = well_style.get('border_color', '#ccc')
                border_style = well_style.get('border_style', '#solid')
                border_size = well_style.get('border_size', 1)
                background = well_style.get('background', color)
                html += f"""
                        <td style="background: {background}; border: {border_size}px {border_style} {border_color}; 
                                 padding: 8px; text-align: center; font-size: 10px; 
                                 width: {box_size}px; height: {int(box_size * 0.8)}px; vertical-align: middle;">
                            {display_value}
                        </td>
                """
            
            html += "</tr>"
        
        html += """
                </tbody>
            </table>
            
            <div style="margin-top: 10px; text-align: center; font-size: 12px; color: #666;">
                Total Wells: 96 | Unique {}: {}
            </div>
        </div>
        """.format(title_suffix, len(unique_values))
        
        return html
    
    def show_medium_layout(self, box_size=50):
        """Display plate colored by medium types"""
        from IPython.display import HTML, display
        display(HTML(self._repr_html_(color_by='medium', box_size=box_size)))
        
    def show_culture_layout(self, box_size=50):
        """Display plate colored by culture types"""
        from IPython.display import HTML, display
        display(HTML(self._repr_html_(color_by='culture', box_size=box_size)))


class PlateSet:

    def __init__(self):
        self.plates = []

    def __iter__(self):
        """Iterate over all plates in the set"""
        for plate in self.plates:
            yield plate

    def get_culture_list(self):
        """Get a list of all unique cultures in the set"""
        return list(set([well.culture for plate in self.plates for well in plate]))
    
    def get_medium_list(self):
        """Get a list of all unique media types in the set"""
        return list(set([well.medium for plate in self.plates for well in plate]))

    def _generate_colors(self, unique_values):
        """Generate distinct colors for different values"""
        colors = [
            '#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7',
            '#DDA0DD', '#98D8C8', '#F7DC6F', '#BB8FCE', '#85C1E9',
            '#F8C471', '#82E0AA', '#F1948A', '#85C1E9', '#D7DBDD',
            '#AED6F1', '#A9DFBF', '#F9E79F', '#D2B4DE', '#AEB6BF'
        ]
        
        color_map = {}
        for i, value in enumerate(unique_values):
            if value == '?':  # Default/empty wells
                color_map[value] = '#F8F9FA'  # Light gray
            else:
                color_map[value] = colors[i % len(colors)]
        
        return color_map

    def add_plate(self, plate):
        """Add a plate to the set"""
        if isinstance(plate, WellPlate96):
            self.plates.append(plate)
        else:
            raise ValueError("Only WellPlate96 objects can be added to PlateSet")

    def __repr_html__(self, color_by='medium', box_size=30, cols=3):
        """
        Generate HTML representation of all plates in a grid layout with unified coloring
        
        Args:
            color_by (str): 'medium' or 'culture' - what to color the wells by
            box_size (int): Size of each well box in pixels (default: 30)
            cols (int): Number of plates per row in the grid (default: 3)
        """
        if not self.plates:
            return "<div style='text-align: center; color: #666;'>No plates in set</div>"
        
        # Collect unique values for coloring across ALL plates
        if color_by == 'medium':
            values = [str(x) for x in self.get_medium_list()]
            title_suffix = "Medium Types"
        elif color_by == 'culture':
            values = [str(x) for x in self.get_culture_list()]
            title_suffix = "Culture Types"
        else:
            raise ValueError("color_by must be 'medium' or 'culture'")
            
        unique_values = sorted(list(set(values)))
        color_map = self._generate_colors(unique_values)
        
        # Start building HTML with unified legend
        html = f"""
        <div style="font-family: Arial, sans-serif; margin: 20px;">
            <h2 style="text-align: center; color: #333;">
                Plate Set - {len(self.plates)} Plates - Colored by {title_suffix}
            </h2>
            
            <!-- Unified Legend -->
            <div style="margin-bottom: 20px; text-align: center;">
                <strong>Legend:</strong><br>
                <div style="display: inline-block; text-align: left;">
        """
        
        # Add legend items
        for value in unique_values:
            color = color_map[value]
            display_name = value if value != '?' else 'Empty/Unknown'
            html += f"""
                    <span style="display: inline-block; margin: 2px 10px;">
                        <span style="width: 15px; height: 15px; background-color: {color}; 
                                   display: inline-block; border: 1px solid #ccc; margin-right: 5px;"></span>
                        {display_name}
                    </span>
            """
        
        html += """
                </div>
            </div>
            
            <!-- Plates Grid -->
            <div style="display: grid; grid-template-columns: repeat({}, 1fr); gap: 20px; justify-items: center;">
        """.format(cols)
        
        # Generate each plate in the grid
        for plate in self.plates:
            html += f"""
                <div style="border: 2px solid #ddd; border-radius: 10px; padding: 15px; background-color: #fafafa;">
                    <h4 style="text-align: center; margin: 0 0 15px 0; color: #555;">
                        {plate.plate_name}
                    </h4>
                    
                    <!-- Plate Layout -->
                    <table style="border-collapse: collapse; margin: 0 auto; border: 2px solid #333;">
                        <thead>
                            <tr>
                                <th style="background-color: #f0f0f0; border: 1px solid #ccc; padding: 4px; font-weight: bold; font-size: 10px;"></th>
            """
            
            # Column headers (1-12)
            for col in plate.cols:
                html += f'<th style="background-color: #f0f0f0; border: 1px solid #ccc; padding: 4px; font-weight: bold; text-align: center; font-size: 10px;">{col}</th>'
            
            html += """
                            </tr>
                        </thead>
                        <tbody>
            """
            
            # Generate plate rows
            for row in plate.rows:
                html += f"""
                            <tr>
                                <td style="background-color: #f0f0f0; border: 1px solid #ccc; padding: 4px; font-weight: bold; text-align: center; font-size: 10px;">{row}</td>
                """
                
                for col in plate.cols:
                    well_id = f"{row}{col}"
                    well = plate.wells[well_id]
                    well_style = plate.wells_style[well_id]

                    if color_by == 'medium':
                        value = str(well.medium)
                    else:
                        value = str(well.culture)
                    
                    color = color_map[value]
                    display_value = value if value != '?' else ''
                    
                    # Truncate long values for display
                    if len(display_value) > 6:
                        display_value = display_value[:4] + '..'

                    border_color = well_style.get('border_color', '#ccc')
                    border_style = well_style.get('border_style', '#solid')
                    border_size = well_style.get('border_size', 1)
                    background = well_style.get('background', color)
                    html += f"""
                            <td style="background: {background}; border: {border_size}px {border_style} {border_color}; 
                                     padding: 2px; text-align: center; font-size: 8px; 
                                     width: {box_size}px; height: {int(box_size * 0.8)}px; vertical-align: middle;">
                                {display_value}
                            </td>
                    """
                
                html += "</tr>"
            
            html += """
                        </tbody>
                    </table>
                </div>
            """
        
        html += """
            </div>
            <div style="margin-top: 20px; text-align: center; font-size: 14px; color: #666;">
                Total Plates: {} | Grid Layout: {} columns | Unique {}: {}
            </div>
        </div>
        """.format(len(self.plates), cols, title_suffix, len(unique_values))
        
        return html

    def show_grid(self, color_by='medium', box_size=50, cols=3):
        """Display all plates in a grid layout"""
        from IPython.display import HTML, display
        display(HTML(self.__repr_html__(color_by=color_by, box_size=box_size, cols=cols)))
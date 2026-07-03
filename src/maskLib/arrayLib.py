# -*- coding: utf-8 -*-
"""
Created 2026

@author: Agrim

arrayLib: parameter-sweep ("dose array") machinery for field-grid layouts.

A chip is treated as a grid of identical fields (e.g. 40 x 40), tiled with
identical tile_nx x tile_ny blocks ("tiles"). A Sweep3D varies parameters
along three axes:

    col  -- varies along x inside a tile (tile_nx steps)
    row  -- varies along y inside a tile (tile_ny steps)
    tile -- varies from tile to tile     (one step per tile)

Each axis is a dict {parameter_name: (start, step)}. The number of steps is
set by the axis (tile_nx, tile_ny, or the number of tiles) and the final
value is computed automatically: final = start + (steps-1)*step. Several
parameters can share one axis (they sweep in lockstep).

Parameters come in two kinds, declared by the registries passed to Sweep3D:
    geometry_params {name: geometry_kwarg} -- change the drawn shapes
    dose_params     {name: LAYER_FAMILY}   -- redirect shapes onto per-dose
                                              layers like BRIDGE_400
Dose layers are generated automatically via Sweep3D.dose_layers().

Every field gets a label like 'A0103' = tile A, column 01, row 03.
Sweep3D.export_workbook() writes an .xlsx with the full parameter table, a
label minimap, and one gradient-colored value map per swept parameter.
"""

import numpy as np


def fmt_value(v):
    '''compact number formatting for layer names and labels (400.0 -> 400)'''
    return '%g' % v


def index_letter(i):
    '''0->A ... 25->Z, 26->AA, ... (spreadsheet style)'''
    s = ''
    while True:
        s = chr(ord('A') + i % 26) + s
        i = i // 26 - 1
        if i < 0:
            return s


def dose_layer(base, params, dose_name):
    '''layer for a structure: BASE if its dose is not swept, else BASE_<dose>'''
    v = params.get(dose_name)
    return base if v is None else '%s_%s' % (base, fmt_value(v))


class Sweep3D:
    '''
    3D parameter sweep over a tiled field grid. See the module docstring for
    the axis conventions. grid_nx/grid_ny are the reference grid (the main
    chiplet); secondary grids (e.g. corner chips) can be queried through
    field(..., grid_nx=..., grid_ny=..., strict=False).
    '''

    def __init__(self, grid_nx, grid_ny, tile_nx, tile_ny,
                 col=None, row=None, tile=None,
                 geometry_params=None, dose_params=None):
        assert grid_nx % tile_nx == 0 and grid_ny % tile_ny == 0, (
            'Tile size %dx%d does not evenly fill the %dx%d field grid -- '
            'pick tile dimensions that divide the grid'
            % (tile_nx, tile_ny, grid_nx, grid_ny))
        self.grid_nx, self.grid_ny = grid_nx, grid_ny
        self.tile_nx, self.tile_ny = tile_nx, tile_ny
        self.tiles_x = grid_nx // tile_nx
        self.tiles_y = grid_ny // tile_ny
        self.n_tiles = self.tiles_x * self.tiles_y
        self.geometry_params = geometry_params or {}
        self.dose_params = dose_params or {}

        # expand each axis {param: (start, step)} into {param: value array};
        # values are rounded so steps like 0.01 stay decimal-exact
        self.axes = {}
        for axname, spec, steps in (('col', col, tile_nx),
                                    ('row', row, tile_ny),
                                    ('tile', tile, self.n_tiles)):
            values = {}
            for pname, (start, step) in (spec or {}).items():
                assert pname in self.geometry_params or pname in self.dose_params, (
                    'unknown sweep parameter: %s (add it to geometry_params '
                    'or dose_params)' % pname)
                values[pname] = np.round(start + step * np.arange(steps), 9)
            self.axes[axname] = values

    # ---- lookups ------------------------------------------------------

    def field(self, ix, iy, grid_nx=None, grid_ny=None, strict=True):
        '''
        (params, label) for field (ix, iy). params is {parameter: value};
        label is e.g. 'A0103' = tile A, column 01, row 03. With strict=True
        the tile must fill the grid evenly; use strict=False for secondary
        chips, where tile parameter values wrap around.
        '''
        gx = self.grid_nx if grid_nx is None else grid_nx
        gy = self.grid_ny if grid_ny is None else grid_ny
        if strict:
            assert gx % self.tile_nx == 0 and gy % self.tile_ny == 0, (
                'Tile size %dx%d does not evenly fill the %dx%d field grid'
                % (self.tile_nx, self.tile_ny, gx, gy))
        ti, tj = ix % self.tile_nx, iy % self.tile_ny
        tile = (ix // self.tile_nx) + (iy // self.tile_ny) * (gx // self.tile_nx)
        params = {}
        for pname, vals in self.axes['col'].items():
            params[pname] = vals[ti]
        for pname, vals in self.axes['row'].items():
            params[pname] = vals[tj]
        for pname, vals in self.axes['tile'].items():
            params[pname] = vals[tile % len(vals)]
        return params, '%s%02d%02d' % (index_letter(tile), ti, tj)

    def param_names(self):
        return sorted({p for ax in self.axes.values() for p in ax})

    def kind(self, pname):
        return 'geometry' if pname in self.geometry_params else 'dose'

    def dose_layers(self):
        '''every dose layer implied by the sweep, e.g. [BRIDGE_400, BRIDGE_440, ...]'''
        out = []
        for axname in ('col', 'row', 'tile'):
            for pname, vals in self.axes[axname].items():
                if pname in self.dose_params:
                    for v in vals:
                        lyr = '%s_%s' % (self.dose_params[pname], fmt_value(v))
                        if lyr not in out:
                            out.append(lyr)
        return out

    # ---- reporting ----------------------------------------------------

    def axis_kind(self, axname):
        '''"Size" (geometry), "Dose", "Mixed", or "None" for one axis'''
        kinds = {('Size' if p in self.geometry_params else 'Dose')
                 for p in self.axes[axname]}
        if not kinds:
            return 'None'
        return kinds.pop() if len(kinds) == 1 else 'Mixed'

    def sweep_type(self):
        '''e.g. "SizeDoseDose" = col sweeps geometry, row and tile sweep dose'''
        return ''.join(self.axis_kind(ax) for ax in ('col', 'row', 'tile'))

    def print_summary(self):
        print('=' * 68)
        print('3D parameter sweep: %s  (col / row / tile)' % self.sweep_type())
        print('Field grid %d x %d; tile %d x %d -> %d x %d = %d tiles (%s-%s)'
              % (self.grid_nx, self.grid_ny, self.tile_nx, self.tile_ny,
                 self.tiles_x, self.tiles_y, self.n_tiles,
                 index_letter(0), index_letter(self.n_tiles - 1)))
        for axname, along, steps in (('col', 'along x in tile', self.tile_nx),
                                     ('row', 'along y in tile', self.tile_ny),
                                     ('tile', 'per tile', self.n_tiles)):
            if not self.axes[axname]:
                continue
            print('%s axis (%d steps, %s):' % (axname.upper(), steps, along))
            for pname, vals in self.axes[axname].items():
                step = vals[1] - vals[0] if len(vals) > 1 else 0
                if self.kind(pname) == 'dose':
                    target = 'dose -> %s_* layers' % self.dose_params[pname]
                else:
                    target = 'geometry -> %s' % self.geometry_params[pname]
                print('    %-20s %8s -> %-8s step %-6s (%s)'
                      % (pname, fmt_value(vals[0]), fmt_value(vals[-1]),
                         fmt_value(step), target))
        n = len(self.dose_layers())
        if n:
            print('Auto-generated dose layers: %d' % n)
        print('=' * 68)

    def legend_lines(self):
        '''compact per-parameter summary for drawing in the chip margin'''
        lines = []
        for tag, axname in (('COLS', 'col'), ('ROWS', 'row'), ('TILES', 'tile')):
            for pname, vals in self.axes[axname].items():
                step = vals[1] - vals[0] if len(vals) > 1 else 0
                lines.append('%s: %s = %s to %s step %s'
                             % (tag, pname, fmt_value(vals[0]),
                                fmt_value(vals[-1]), fmt_value(step)))
        return lines

    # ---- export -------------------------------------------------------

    def export_workbook(self, path):
        '''
        Write an .xlsx workbook (requires openpyxl):
          'parameters'  -- label / ix / iy / every parameter value per field
          'map'         -- field labels laid out like the chip
                           (top spreadsheet row = top of chip)
          one sheet per swept parameter -- its values laid out like the
          chip, with a min->max color gradient so the sweep is visible
        '''
        from openpyxl import Workbook
        from openpyxl.formatting.rule import ColorScaleRule
        from openpyxl.utils import get_column_letter

        pnames = self.param_names()
        wb = Workbook()
        ws = wb.active
        ws.title = 'parameters'
        ws.append(['label', 'ix', 'iy'] + pnames)
        for ix in range(self.grid_nx):
            for iy in range(self.grid_ny):
                params, flabel = self.field(ix, iy)
                ws.append([flabel, ix, iy] + [float(params[p]) for p in pnames])

        def grid_sheet(title, cellvalue):
            s = wb.create_sheet(title)
            s.append(['iy\\ix'] + list(range(self.grid_nx)))
            for iy in reversed(range(self.grid_ny)):
                s.append([iy] + [cellvalue(ix, iy) for ix in range(self.grid_nx)])
            return s

        grid_sheet('map', lambda ix, iy: self.field(ix, iy)[1])

        # one gradient-colored value map per swept parameter (sheet titles
        # are capped at Excel's 31-character limit)
        data_range = 'B2:%s%d' % (get_column_letter(self.grid_nx + 1),
                                  self.grid_ny + 1)
        for pname in pnames:
            s = grid_sheet(pname[:31],
                           lambda ix, iy, p=pname: float(self.field(ix, iy)[0][p]))
            s.conditional_formatting.add(data_range, ColorScaleRule(
                start_type='min', start_color='FFFFFFFF',
                end_type='max', end_color='FFFF4444'))
        wb.save(path)
        return path

from __future__ import annotations

import math
import typing as t

import numpy as np


def extract_wulff_geometry(wulff: t.Any, rounding: int = 8) -> dict[str, t.Any]:
    """
    Convert a Pymatgen WulffShape into serializable geometry.

    Args:
        wulff: Pymatgen WulffShape instance.
        rounding: Decimal places used to deduplicate vertices.

    Returns:
        Dict with vertices (list of [x, y, z]), faces (list of index lists),
        facet_normals, all_points (raw wulff_pt_list), and bounding_radius.
    """
    vertices: list[list[float]] = []
    vertex_index: dict[tuple[float, ...], int] = {}
    faces: list[list[int]] = []
    facet_normals: list[tuple[float, float, float]] = []

    for facet in getattr(wulff, 'facets', []):
        points = getattr(facet, 'points', None)
        if not points:
            continue

        normal = getattr(facet, 'normal', None)
        if normal is not None:
            facet_normals.append(tuple(float(x) for x in normal))

        for polygon in points:
            if not polygon:
                continue

            face: list[int] = []
            for point in polygon:
                key = tuple(round(float(coord), rounding) for coord in point)
                idx = vertex_index.get(key)
                if idx is None:
                    idx = len(vertices)
                    vertex_index[key] = idx
                    vertices.append([float(c) for c in point])
                face.append(idx)

            if len(face) >= 3:
                faces.append(face)

    all_points = [
        [float(coord) for coord in point]
        for point in getattr(wulff, 'wulff_pt_list', [])
    ]
    bounding_radius = max(
        (math.sqrt(np.dot(point, point)) for point in all_points), default=0.0
    )

    return {
        'vertices': vertices,
        'faces': faces,
        'facet_normals': facet_normals,
        'points': all_points,
        'bounding_radius': bounding_radius,
    }


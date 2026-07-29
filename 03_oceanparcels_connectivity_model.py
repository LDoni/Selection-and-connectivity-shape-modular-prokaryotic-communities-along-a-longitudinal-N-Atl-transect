from __future__ import annotations

import argparse
import json
import shutil
from datetime import timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

import parcels
from parcels.tools.statuscodes import StatusCode


# Follows the official OceanParcels quickstart pattern:
# https://docs.oceanparcels.org/en/stable/examples/parcels_tutorial.html
# FieldSet -> ParticleSet.from_list -> ParticleFile -> pset.execute(AdvectionRK4).
#
# Output is read as Zarr/xarray following:
# https://docs.oceanparcels.org/en/stable/examples/tutorial_output.html

EARTH_RADIUS_KM = 6371.0088


def read_stations(path: str) -> pd.DataFrame:
    stations = pd.read_csv(path).rename(
        columns={
            "Nominativo.campione.": "station",
            "data.x": "date",
            "Latitude.x": "lat",
            "Longitude.x": "lon",
            "sample_id": "station",
            "sampling_date": "date",
            "latitude": "lat",
            "longitude": "lon",
        }
    )
    stations["date"] = pd.to_datetime(stations["date"], dayfirst=True)
    return stations.loc[:, ["station", "date", "lat", "lon"]].reset_index(drop=True)


def particle_cloud(lon: float, lat: float, side: int, jitter_deg: float):
    offsets = np.linspace(-jitter_deg, jitter_deg, side)
    return [(lon + dx, lat + dy) for dx in offsets for dy in offsets]


def make_fieldset(currents: str):
    filenames = {"U": currents, "V": currents}
    variables = {"U": "uo", "V": "vo"}
    dimensions = {
        "U": {"lon": "longitude", "lat": "latitude", "depth": "depth", "time": "time"},
        "V": {"lon": "longitude", "lat": "latitude", "depth": "depth", "time": "time"},
    }
    return parcels.FieldSet.from_netcdf(filenames, variables, dimensions, mesh="spherical")


def run_station(fieldset, station: pd.Series, direction: str, outdir: Path, days: int, dt_hours: float, output_hours: float, side: int, jitter_deg: float):
    cloud = particle_cloud(station.lon, station.lat, side, jitter_deg)
    lon = np.array([p[0] for p in cloud], dtype=np.float32)
    lat = np.array([p[1] for p in cloud], dtype=np.float32)
    time = np.array([station.date.to_datetime64()] * len(cloud), dtype="datetime64[ns]")

    pset = parcels.ParticleSet.from_list(
        fieldset=fieldset,
        pclass=parcels.JITParticle,
        lon=lon,
        lat=lat,
        time=time,
    )

    zarr_path = outdir / "trajectories" / f"{direction}_{station.station}.zarr"
    if zarr_path.exists():
        shutil.rmtree(zarr_path)

    output_file = pset.ParticleFile(
        name=str(zarr_path),
        outputdt=timedelta(hours=output_hours),
        chunks=(len(cloud), 256),
    )
    output_file.metadata["source_station"] = str(station.station)
    output_file.metadata["direction"] = direction

    dt = timedelta(hours=dt_hours)
    if direction == "backward":
        dt = -dt

    def DeleteOutOfBounds(particle, fieldset, time):
        if particle.state >= StatusCode.ErrorOutOfBounds:
            particle.delete()

    pset.execute(
        parcels.AdvectionRK4 + pset.Kernel(DeleteOutOfBounds),
        runtime=timedelta(days=days),
        dt=dt,
        output_file=output_file,
    )
    return zarr_path


def haversine_km(lat1, lon1, lat2, lon2):
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return EARTH_RADIUS_KM * 2 * np.arcsin(np.sqrt(a))


def summarize_station_zarr(zarr_path: Path, source_index: int, stations: pd.DataFrame, release_time, hit_radius_km: float):
    ds = xr.open_zarr(zarr_path)
    lon = np.asarray(ds["lon"], dtype=float)
    lat = np.asarray(ds["lat"], dtype=float)
    time = np.asarray(ds["time"])
    release_time = np.datetime64(release_time)
    age_days = np.abs((time - release_time) / np.timedelta64(1, "D"))

    conn = np.zeros(len(stations), dtype=float)
    min_days = np.full(len(stations), np.nan)

    for sink_index, sink in stations.iterrows():
        d = haversine_km(lat, lon, sink.lat, sink.lon)
        hit = d <= hit_radius_km
        conn[sink_index] = hit.any(axis=1).mean()
        if hit.any():
            min_days[sink_index] = np.nanmin(age_days[hit])

    ds.close()
    return conn, min_days


def build_matrices(outdir: Path, stations: pd.DataFrame, direction: str, hit_radius_km: float):
    names = stations["station"].to_list()
    conn = np.zeros((len(stations), len(stations)), dtype=float)
    min_days = np.full((len(stations), len(stations)), np.nan)

    for source_index, station in stations.iterrows():
        zarr_path = outdir / "trajectories" / f"{direction}_{station.station}.zarr"
        conn[source_index, :], min_days[source_index, :] = summarize_station_zarr(
            zarr_path,
            source_index,
            stations,
            station.date.to_datetime64(),
            hit_radius_km,
        )

    return (
        pd.DataFrame(conn, index=names, columns=names),
        pd.DataFrame(min_days, index=names, columns=names),
    )


def matrix_to_edges(conn: pd.DataFrame, min_days: pd.DataFrame, stations: pd.DataFrame):
    coord = stations.set_index("station")
    rows = []
    for source in conn.index:
        for sink in conn.columns:
            if source == sink or conn.loc[source, sink] <= 0:
                continue
            rows.append(
                {
                    "source": source,
                    "sink": sink,
                    "connectivity": conn.loc[source, sink],
                    "min_days": min_days.loc[source, sink],
                    "source_lon": coord.loc[source, "lon"],
                    "source_lat": coord.loc[source, "lat"],
                    "sink_lon": coord.loc[sink, "lon"],
                    "sink_lat": coord.loc[sink, "lat"],
                }
            )
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stations", required=True)
    parser.add_argument("--currents", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--days", type=int, default=365)
    parser.add_argument("--dt-hours", type=float, default=6)
    parser.add_argument("--output-hours", type=float, default=6)
    parser.add_argument("--cloud-side", type=int, default=7)
    parser.add_argument("--jitter-deg", type=float, default=0.08)
    parser.add_argument("--hit-radius-km", type=float, default=50)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    (outdir / "trajectories").mkdir(parents=True, exist_ok=True)

    stations = read_stations(args.stations)
    fieldset = make_fieldset(args.currents)

    for direction in ["forward", "backward"]:
        for _, station in stations.iterrows():
            print(f"running {direction} {station.station}", flush=True)
            run_station(
                fieldset,
                station,
                direction,
                outdir,
                args.days,
                args.dt_hours,
                args.output_hours,
                args.cloud_side,
                args.jitter_deg,
            )

    forward_conn, forward_time = build_matrices(outdir, stations, "forward", args.hit_radius_km)
    backward_conn, backward_time = build_matrices(outdir, stations, "backward", args.hit_radius_km)

    forward_conn.to_csv(outdir / "forward_connectivity_matrix.csv")
    forward_time.to_csv(outdir / "forward_min_travel_days.csv")
    backward_conn.to_csv(outdir / "backward_connectivity_matrix.csv")
    backward_time.to_csv(outdir / "backward_min_travel_days.csv")
    matrix_to_edges(forward_conn, forward_time, stations).to_csv(outdir / "forward_connected_edges.csv", index=False)
    matrix_to_edges(backward_conn, backward_time, stations).to_csv(outdir / "backward_connected_edges.csv", index=False)

    metadata = {
        "software": "OceanParcels",
        "parcels_version": parcels.__version__,
        "fieldset": "parcels.FieldSet.from_netcdf",
        "particleset": "parcels.ParticleSet.from_list",
        "kernel": "parcels.AdvectionRK4",
        "output": "pset.ParticleFile zarr; read with xarray.open_zarr",
        "days_per_station": args.days,
        "dt_hours": args.dt_hours,
        "output_hours": args.output_hours,
        "cloud_side": args.cloud_side,
        "particles_per_station": args.cloud_side**2,
        "n_stations": int(len(stations)),
        "release_events_per_station_and_direction": 1,
        "release_schedule": "one particle cloud at each station sampling date",
        "total_particles_per_direction": int(len(stations) * args.cloud_side**2),
        "jitter_deg": args.jitter_deg,
        "hit_radius_km": args.hit_radius_km,
        "connectivity_definition": "fraction of the source cloud entering the sink radius at least once",
        "connectivity_denominator": args.cloud_side**2,
        "multiple_hit_handling": "one particle counts once per sink; it may count once for each distinct sink reached",
    }
    (outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
    
    
    
    
    
    
 
    
    

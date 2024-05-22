import dask.array as da
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader
from ome_zarr.writer import write_image
from ome_zarr.scale import Scaler
import fsspec
import zarr
from dask.diagnostics import ProgressBar

in_zarr_url = snakemake.params.in_zarr
loc = parse_url(in_zarr_url,mode='r')
zarr_reader = Reader(loc).zarr


attrs=zarr_reader.root_attrs


darr = zarr_reader.load(0)[:,
                snakemake.params.fov['z'][0]:snakemake.params.fov['z'][1],
                snakemake.params.fov['y'][0]:snakemake.params.fov['y'][1],
                snakemake.params.fov['x'][0]:snakemake.params.fov['x'][1]]


print(f'darr shape: {darr.shape}')
group = zarr.open_group(snakemake.output.ome_zarr,path='/',mode='w')
max_layer =5
scaling_method = 'local_mean'
scaler = Scaler(max_layer=max_layer,method=scaling_method)
print(attrs)

coordinate_transformations = [dataset['coordinateTransformations'] for dataset in attrs['multiscales'][0]['datasets'] ]
omero = attrs['omero']
axes = attrs['multiscales'][0]['axes']

with ProgressBar():
    write_image(image=darr.rechunk(chunks=(1,100,200,200)),
                            group=group,
                            scaler=scaler,
                            coordinate_transformations=coordinate_transformations,
                            axes=axes,
                            storage_options={'dimension_separator': '/'},
                            metadata={'omero':omero}
                                )



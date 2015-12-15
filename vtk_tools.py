import vtk as v


def save_lesion(output_file, input_file, field, limits):
    reader = v.vtkXMLUnstructuredGridReader()
    reader.SetFileName(input_file)
    reader.Update()

    threshold = v.vtkThreshold()
    threshold.SetInput(reader.GetOutput())

    if limits[0] is None:
        threshold.ThresholdByLower(limits[1])
    elif limits[1] is None:
        threshold.ThresholdByUpper(limits[0])
    else:
        threshold.ThresholdBetween(*limits)

    threshold.SetInputArrayToProcess(0, 0, 0, v.vtkDataObject.FIELD_ASSOCIATION_CELLS, field)
    threshold.Update()

    extract_surface = v.vtkDataSetSurfaceFilter()
    extract_surface.SetInput(threshold.GetOutput())
    extract_surface.Update()

    writer = v.vtkXMLPolyDataWriter()
    writer.SetFileName(output_file)
    writer.SetInput(extract_surface.GetOutput())
    writer.Write()

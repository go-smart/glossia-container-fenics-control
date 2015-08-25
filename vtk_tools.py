import vtk as v


def save_lesion(output_file, input_file, field, limits):
    reader = v.vtkXMLUnstructuredGridReader()
    reader.SetFileName(input_file)
    reader.Update()

    threshold = v.vtkThreshold()
    threshold.SetInputData(reader.GetOutput())
    threshold.SetInputArrayToProcess(0, 0, 0, field, threshold.vtkDataSetAttributes.SCALARS)

    if limits[0] is None:
        threshold.ThresholdByLower(limits[1])
    elif limits[1] is None:
        threshold.ThresholdByUpper(limits[0])
    else:
        threshold.ThresholdBetween(*limits)

    extract_surface = v.vtkDataSetSurfaceFilter()
    extract_surface.SetInputData(threshold.GetOutput())

    writer = v.vtkXMLPolyDataWriter()
    writer.SetFileName(output_file)
    writer.SetInputData(extract_surface.GetOutput())
    writer.Write()

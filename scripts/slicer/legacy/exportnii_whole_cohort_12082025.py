outputFolder = "/media/jezdip1/UBURYZEN24_48T1/motol/AZV_PET/normy/whole_cohort_12082025/"
print(f'reading database')
import csv
import traceback

from DICOMLib import DICOMUtils
patientUIDs = slicer.dicomDatabase.patients()

#with open('/home/jezdip1/motol/ISARGUBU/media/shared_storage/motol/AZV_PET/normy/PET_Normy_09102023.csv', mode='r') as file:
with open('/media/jezdip1/UBURYZEN24_48T1/motol/AZV_PET/normy/all_series_metadata_export_12082025_pje.csv', mode='r') as file:
    reader = csv.reader(file)
    # Přeskočení hlavičky, pokud je v souboru
    next(reader, None)
    # Procházení každého řádku v CSV souboru
    for row in reader:
        # Získání hodnot z sloupců 'ID' a 'UNIS'
        ID = row[3]
        UNIS=row[2]
        # Zpracování dat, např. výpis ID a UNIS
        print("ID:", ID)
        print("UNIS:", UNIS)
        print("-----")  # Oddělovač mezi záznamy
        try:
            print(f'trying to load patient by patient id ',ID)
            loadedNodeIDs = DICOMUtils.loadPatientByPatientID(ID)
            print(f'loaded patient id ',ID)
            for loadedNodeID in loadedNodeIDs:
        	# Check if we want to save this node
                node = slicer.mrmlScene.GetNodeByID(loadedNodeID)
        # Only export images
                if not node or not node.IsA('vtkMRMLScalarVolumeNode'):
                    continue
        # Construct filename
                shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
                seriesItem = shNode.GetItemByDataNode(node)
                studyItem = shNode.GetItemParent(seriesItem)
                patientItem = shNode.GetItemParent(studyItem)
#            filename = shNode.GetItemAttribute(patientItem, 'DICOM.PatientID')
                filename = UNIS
                filename += '_' + shNode.GetItemAttribute(studyItem, 'DICOM.StudyDate')
                filename += '_' + shNode.GetItemAttribute(seriesItem, 'DICOM.SeriesNumber')
                filename += '_' + shNode.GetItemAttribute(seriesItem, 'DICOM.Modality')
                filename += '_' + shNode.GetItemAttribute(seriesItem, 'DICOM.SeriesDescription')
                filename = slicer.app.ioManager().forceFileNameValidCharacters(filename) + ".nii"
        # Save node
                print(f'Write {node.GetName()} to {filename}')
                success = slicer.util.saveNode(node, outputFolder+filename)
            slicer.mrmlScene.Clear()
        except:
            print(f"Something wrong with the patient {ID}, most probably does not exists in database")
#            traceback.print_exc()

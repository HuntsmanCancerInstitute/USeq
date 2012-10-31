package edu.utah.bass;

import java.awt.Component;
import java.io.File;
import java.io.IOException;
import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.freehep.util.export.ExportDialog;
import org.freehep.util.export.ExportFileType;

/**
 * Prints a component to a file.
 * That file can be in any of the formats supported by the FreeHep VectorGraphics library: http://java.freehep.org/vectorgraphics
 * These include:  PostScript, PDF, EMF, SVF, Flash SWF, GIF, PNG, JPG and PPM.
 * Turns off/on double-buffering in Swing and NGSDK components.
 *
 */
public final class ComponentWriter {

	/** Show the export dialog that allows exporting in a variety of graphics
	 *  formats using the FreeHep libraries.
	 */
	public static void showExportDialog(Component c) {
		ExportDialog  export = new ExportDialog();
		export.showExportDialog(c, "Export view as ...", c, "export");
	}

	/** Show the export dialog that allows exporting in a variety of graphics
	 *  formats using the FreeHep libraries.
	 */
	public static boolean exportComponent(File f, Component c, ExportFileType eft) {
		ExportNoDialog exportNoDialog = new ExportNoDialog(f);
		return exportNoDialog.exportIt(f, c, eft);
	}

	private final static class ExportNoDialog extends ExportDialog {
		File f = null;
		Properties props = new Properties();
		ExportNoDialog(File f) {
			this.f = f;
		}

		@Override
		protected String selectFile() {
			return f.getAbsolutePath();
		}

		/** Called to actually write out the file.
		 * Override this method to provide special handling
		 * @return true if the file was written, or false to cancel operation
		 */
		@Override
		protected boolean writeFile(Component component, ExportFileType t) throws IOException {
			t.exportToFile(f, component, this, props, "");
			//props.put(SAVE_AS_FILE, file.getText());
			//props.put(SAVE_AS_TYPE, currentType().getFileFilter().getDescription());
			return true;
		}
		boolean exportIt(File f, Component c, ExportFileType eft) {
			try {
				this.selectFile();
				return this.writeFile(c, eft);
			} catch (IOException ex) {
				Logger.getLogger(ComponentWriter.class.getName()).log(Level.SEVERE, null, ex);
				return false;
			}
		}
	}
}


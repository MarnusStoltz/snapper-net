package beast.app.draw;
 
import java.util.List;

import javax.swing.Box;

import beast.app.beauti.BeautiDoc;
import beast.app.draw.ListInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import snappNetProject.core.SnappNetPrior;

public class SNAPPPriorListEditor extends ListInputEditor {
    public SNAPPPriorListEditor(BeautiDoc doc) {
		super(doc);
	}

	private static final long serialVersionUID = 1L;

    public Class<?> baseType() {
        return SnappNetPrior.class;
    }
    
    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpand, boolean bAddButtons) {
		m_bAddButtons = bAddButtons;
    	m_bExpandOption = bExpand;
        m_input = input;
        m_beastObject = plugin;
		this.itemNr = itemNr;

        m_listBox = Box.createVerticalBox();
        // list of inputs 
        for (Object o : (List<?>) input.get()) {
        	
        		System.out.println("CE voici l objet ds snapriorlisteditor"+ o +"\n");
        	
            if (o instanceof SnappNetPrior) {
            	SnappNetPrior plugin2 = (SnappNetPrior) o;
            	doc.getInputEditorFactory().addInputs(m_listBox, plugin2, this, null, doc);
            }
        }
		add(m_listBox);
        updateState();
    }

}

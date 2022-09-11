package beast.app.draw;
 
import java.util.List;

import javax.swing.Box;
import javax.swing.JButton;

import beast.app.beauti.BeautiDoc;
import beast.app.draw.InputEditor;
import beast.app.draw.ListInputEditor;
import beast.app.draw.ParameterInputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.sitemodel.SiteModel;
import snappNetProject.core.SnapData;
import snappNetProject.core.SnappNetNetworkLikelihood;
import snappNetProject.core.SnappNetSubstitutionModel;
 
public class SNAPPTreeLikelihoodEditor extends ListInputEditor {
    public SNAPPTreeLikelihoodEditor(BeautiDoc doc) {
		super(doc);
	}

	private static final long serialVersionUID = 1L;

    public Class<?> baseType() {
        return SnappNetNetworkLikelihood.class;
    }
    
    SnappNetSubstitutionModel substModel;
    SnapData data;
    JButton muButton;
    
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
        		System.out.println("CEE voici l objet ds snaptreelikelihood"+ o +"\n");
            if (o instanceof SnappNetNetworkLikelihood) {
            	SnappNetNetworkLikelihood plugin2 = (SnappNetNetworkLikelihood) o;
            	if (plugin2.dataInput.get() instanceof SnapData) {
            		data = (SnapData) plugin2.dataInput.get();
            	}
            	substModel = (SnappNetSubstitutionModel) ((SiteModel.Base) plugin2.siteModelInput.get()).substModelInput.get();
            	doc.getInputEditorFactory().addInputs(m_listBox, substModel, this, null, doc);
            	doc.getInputEditorFactory().addInputs(m_listBox, plugin2, this, null, doc);
            	muButton = new JButton("Calc mutation rates");
            	muButton.setToolTipText("Calcaulate mutation rates based on data in the alignment");
            	muButton.addActionListener(e -> setUpMutationRates());
            	add(muButton);
            }
        }
		add(m_listBox);
        updateState();
    }
    

    private Object setUpMutationRates() {
    	double proportionZeros = data.getProportionZeros();
    	double muU = 1 / (2.0 * (1.0 - proportionZeros));
    	double muV = 1 / (2.0 * proportionZeros);
    	RealParameter pU = substModel.m_pU.get();
    	pU.valuesInput.setValue(muU + "", pU);
    	RealParameter pV = substModel.m_pV.get();
    	pV.valuesInput.setValue(muV + "", pV);
    	refreshPanel();
		return null;
	}


	public InputEditor createMutationRateVEditor() throws Exception {
    	ParameterInputEditor editor = (ParameterInputEditor) doc.getInputEditorFactory().createInputEditor(substModel.m_pV, substModel, doc);
    	editor.m_isEstimatedBox.setVisible(false);
    	return editor;
    }

}

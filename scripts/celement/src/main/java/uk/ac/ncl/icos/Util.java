package uk.ac.ncl.icos;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.net.URI;

import org.sbolstandard.core.DnaSequence;
import org.sbolstandard.core2.AccessType;
import org.sbolstandard.core2.Component;
import org.sbolstandard.core2.ComponentDefinition;
import org.sbolstandard.core2.DirectionType;
import org.sbolstandard.core2.FunctionalComponent;
import org.sbolstandard.core2.Interaction;
import org.sbolstandard.core2.ModuleDefinition;
import org.sbolstandard.core2.OrientationType;
import org.sbolstandard.core2.Participation;
import org.sbolstandard.core2.SBOLDocument;
import org.sbolstandard.core2.SBOLWriter;
import org.sbolstandard.core2.Sequence;
import org.sbolstandard.core2.SequenceAnnotation;
import org.sbolstandard.core2.SystemsBiologyOntology;


public class Util {

	
	public static void write(SBOLDocument doc, File file) throws Exception
	{
		try
		{
			OutputStream stream=new FileOutputStream(file);
			SBOLWriter.write(doc, stream);
			stream.close();			
		}
		catch (Exception e)
		{
			throw new Exception("Could not read the sbol file",e);
		}
	}
	
	public static String write(SBOLDocument doc) throws Exception
	{
		try
		{
			OutputStream stream=new ByteArrayOutputStream();
			SBOLWriter.write(doc, stream);
			String output=stream.toString();			
			return output;
		}
		catch (Exception e)
		{
			throw new Exception("Could not read the sbol file",e);
		}
	}
	
	 public static Interaction createPromoterInduction(ModuleDefinition moduleDef,ComponentDefinition promoterCompDef, ComponentDefinition tfCompDef)
	    {
	    	return createPromoterInduction(moduleDef, promoterCompDef, tfCompDef, null);
	    }
	    
	    
	    public static Interaction createPromoterInduction(ModuleDefinition moduleDef,ComponentDefinition promoter, ComponentDefinition tf, URI tfRole)
	    {
	    	FunctionalComponent fcPromoter=createFunctionalComponent(moduleDef, promoter);
	    	FunctionalComponent fcTF= createFunctionalComponent(moduleDef, tf);
	    	
	    	String interactionId=tf.getDisplayId() + "activates" + promoter.getDisplayId();
	    	Interaction interaction= moduleDef.getInteraction(interactionId);
	    	
	    	if (interaction==null){
		    	interaction=moduleDef.createInteraction(interactionId, SystemsBiologyOntology.STIMULATION);
		    	
		    	Participation promoterParticipation=interaction.createParticipation("promoter_participation", fcPromoter.getIdentity());
		    	promoterParticipation.addRole(SystemsBiologyOntology.PROMOTER);
		    	
		    	Participation tfParticipation=interaction.createParticipation("tf_participation", fcTF.getIdentity());
		    	tfParticipation.addRole(SystemsBiologyOntology.STIMULATOR);
		    	if (tfRole!=null)
		    	{
		    		tfParticipation.addRole(tfRole);
		    	}
		    }   	
	    	return interaction;
	    }
	    
	    public static Interaction createPromoterRepression(ModuleDefinition moduleDef,ComponentDefinition promoter, ComponentDefinition tf)
	    {
	    	FunctionalComponent fcPromoter=createFunctionalComponent(moduleDef, promoter);
	    	FunctionalComponent fcTF= createFunctionalComponent(moduleDef, tf);
	    	
	    	String interactionId=tf.getDisplayId() + "represses" + promoter.getDisplayId();
	    	Interaction interaction= moduleDef.getInteraction(interactionId);
	    	
	    	if (interaction==null){
		    	interaction=moduleDef.createInteraction(interactionId, SystemsBiologyOntology.INHIBITION);
		    	
		    	Participation promoterParticipation=interaction.createParticipation("promoter_participation", fcPromoter.getIdentity());
		    	promoterParticipation.addRole(SystemsBiologyOntology.PROMOTER);
		    	
		    	Participation tfParticipation=interaction.createParticipation("tf_participation", fcTF.getIdentity());
		    	tfParticipation.addRole(SystemsBiologyOntology.INHIBITOR);
		    	/*if (tfRole!=null)
		    	{
		    		tfParticipation.addRole(tfRole);
		    	}*/
		    }   	
	    	return interaction;
	    }
	    
	    
	    private static FunctionalComponent createFunctionalComponent(ModuleDefinition moduleDef,ComponentDefinition compDef)
	    {
	    	FunctionalComponent fc=moduleDef.getFunctionalComponent(compDef.getDisplayId());
	    	if (fc==null)
	    	{
	    			fc=moduleDef.createFunctionalComponent(compDef.getDisplayId(), AccessType.PUBLIC, compDef.getIdentity(), DirectionType.INOUT);
	    	}
	    	return fc;
	    }
	    
	    public static Interaction createTranslationInteraction(ModuleDefinition moduleDef,ComponentDefinition cds, ComponentDefinition protein) 
	    {	    	
		    	FunctionalComponent fcCds=createFunctionalComponent(moduleDef, cds);
		    	FunctionalComponent fcProtein=createFunctionalComponent(moduleDef, protein);
		    	
		    	Interaction interaction=moduleDef.getInteraction(protein.getDisplayId() + "translation");
		    	if (interaction==null)
		    	{
		    			interaction=moduleDef.createInteraction(protein.getDisplayId() + "translation", SystemsBiologyOntology.TRANSLATION);
		    			Participation cdsParticipation=interaction.createParticipation("cds_participation", fcCds.getIdentity());
		    	    	cdsParticipation.addRole(SystemsBiologyOntology.MODIFIER);
		    	    	
		    	    	Participation proteinParticipation=interaction.createParticipation("tf_participation", fcProtein.getIdentity());
		    	    	proteinParticipation.addRole(SystemsBiologyOntology.PRODUCT);	    	
		    	}
		    	
		    	return interaction;	    	
	    }

		public static ComponentDefinition createDNAComponent(SBOLDocument doc, String id, String name, String description, URI role) {
			ComponentDefinition compDef=doc.createComponentDefinition(id, ComponentDefinition.DNA);
			compDef.setName(name);
			compDef.setDescription(description);			
			compDef.addRole(role);
			return compDef;
		}
		
		public static ComponentDefinition createDNAComponent(SBOLDocument doc, String id, String name, String description, URI role, String sequence) {
			ComponentDefinition compDef=createDNAComponent(doc, id, name, description, role);
			Sequence seq=createDNASequence(doc, compDef, sequence);
			compDef.addSequence(seq);
			return compDef;
		}
		
		public static ComponentDefinition createProteinComponent(SBOLDocument doc, String id, String name, String description, URI role, String sequence) {
			ComponentDefinition compDef=doc.createComponentDefinition(id, ComponentDefinition.PROTEIN);
			compDef.setName(name);
			compDef.setDescription(description);	
			if (role!=null)
			{
				compDef.addRole(role);
			}
			if (sequence!=null)
			{
				Sequence seq=createProteinSequence(doc, compDef, sequence);
				compDef.addSequence(seq);
			}
			return compDef;
		}

		private static Sequence createDNASequence(SBOLDocument doc, ComponentDefinition compDef, String sequence)
		{
			return doc.createSequence(compDef.getDisplayId() + "_seq", sequence,URI.create("http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html"));
		}
		
		private static Sequence createProteinSequence(SBOLDocument doc, ComponentDefinition compDef, String sequence)
		{
			return doc.createSequence(compDef.getDisplayId() + "_seq", sequence,URI.create("http://www.chem.qmul.ac.uk/iubmb/misc/aa.html"));//TODO Use the correct one
		}
		
		public static void add(SBOLDocument doc, ComponentDefinition design, ComponentDefinition subComponentDefinition) {
			Component subComponent=design.createComponent(design.getDisplayId() + "_" + subComponentDefinition.getDisplayId(), AccessType.PUBLIC, subComponentDefinition.getIdentity());
			String subComponentSequence=subComponentDefinition.getSequences().iterator().next().getElements();
			int start=1;
			int end=1;
			if (design.getSequences()==null || design.getSequences().size()==0)			
			{
				design.addSequence(createDNASequence(doc, design , subComponentSequence));
				start=1;
			
			}
			else
			{
				Sequence seq=design.getSequences().iterator().next();
				String designSequence=seq.getElements();
				seq.setElements(designSequence + subComponentSequence);
				start=designSequence.length() +1;
			}
			end=start + subComponentSequence.length() -1;			
			SequenceAnnotation anno=design.createSequenceAnnotation("anno" + start + "_" + end, "location" + start  + "_" + end,start, end, OrientationType.INLINE);
			anno.setComponent(subComponent.getIdentity());
		}
		
	    
}

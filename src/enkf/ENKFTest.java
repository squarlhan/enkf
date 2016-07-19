package enkf;

import java.io.File;

import org.openda.interfaces.IStochModelFactory;
import org.openda.interfaces.IStochObserver;
import org.openda.interfaces.IVector;
import org.openda.models.oscillator.OscillatorStochModelFactory;
import org.openda.utils.CsvStochObserver;
import org.openda.utils.StochVector;

public class ENKFTest {
	
	private long initSeed=1234567890;
	private File testRunDataDir = new File("C:/Users/install/Desktop/hxs/workspace/openda_public/opendaTestRuns/algorithms/org/openda/algorithms");
	
	public void testEnKF() {
		System.out.println("========================================================");
		System.out.println(" Test EnKF (Classical Ensemble Kalman Filter.");
		System.out.println("========================================================");
		//generate observations
		StochVector.setSeed(initSeed);
		IStochModelFactory factory = new OscillatorStochModelFactory();
		factory.initialize(null, new String[]{"<oscillatorConfig><simulationTimespan>[0.0,0.05,10.0]</simulationTimespan><parameters names=\"x1,x2\">[4,2]</parameters><parameterUncertainty names=\"x1,x2\">[0.1,0.1]</parameterUncertainty><systemNoise>{[0.0,0.0],[0.3,0.3]}</systemNoise><initialState>[4,2]</initialState><initialStateUncertainty>[0.8,0.8]</initialStateUncertainty></oscillatorConfig>"});
		IStochObserver obsGenerated = new CsvStochObserver();
		obsGenerated.initialize(testRunDataDir, new String[]{"observations_oscillator_generated.csv"});
		//create and run filter
		AbstractSequentialEnsembleAlgorithm enkf = new EnKF();
		enkf.initialize(null, new String[]{"<algorithm><ensembleModel stochParameter=\"false\" stochForcing=\"true\" stochInit=\"true\" /></algorithm>"});
		enkf.setStochComponents(obsGenerated,factory);
		enkf.prepare();
		enkf.run();
		// state at final time
		IVector x = ((EnKF) enkf).getCurrentState();
		System.out.println("x = "+x);
		System.out.println("Should be close to x = [-0.2315739613743379,0.028205645374715847]");

	}
	
	public static void main(String[] args){
		
		ENKFTest kft = new ENKFTest();
		kft.testEnKF();
	}

}

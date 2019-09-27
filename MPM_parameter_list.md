# MPM sequence 

## Specific  acquisition parameters

This list the parameters as being used by the [hMRI toolbox]([http://hmri.info](http://www.hmri.info/)), as of September 2019.

---


| Fieldnames queried by the hMRI toolbox |Fieldnames for BIDS | Already existing in BIDS? | Description/comments |
| ---- | ---- | ---- | ---- |
| **GENERAL** | | | |
| RepetitionTime | RepetitionTimeExcitation	| y |  |
| EchoTime | EchoTime | y |  |
| FlipAngle | FlipAngle | y |  |
| **NAMING** | | | |
| ScanningSequence | ScanningSequence | y | |
| SequenceName | SequenceName | y | |
| ProtocolName | ProtocolName | N | What does it correspond to exactly? Add def in BIDS? |
| ManufacturerModelName | ManufacturerModelName | y |  |
| **RF SPOILING** | | | |
| RFSpoilingPhaseIncrement | SpoilingRFPhaseIncrement | y |  |
| spoilingGradientMoment | SpoilingGradientMoment | y | |
| spoilingGradientDuration | SpoilingGradientDuration | y | |
| **B1 MAPPING** | | ||
| B1mapNominalFAValues | B1mapNominalFAValues | N | Extend `FlipAngle` with `NominalFlipAngle`? |
| B1mapMixingTime | MixingTime | N | Anything close? |
| epiReadoutDuration | epiReadoutTime | N | Is this equivalent to `TotalReadoutTime`? |
| PhaseEncodingDirectionSign | PhaseEncodingDirectionSign | N | Same as `SliceEncodingDirection`? |
| BandwidthPerPixelRO | BandwidthPerPixelRO | N | Seems linked to `EffectiveEchoSpacing` but? |
| PELinesPAT | PELinesPAT | N | ?? |
| NumberOfMeasurements | NumberOfMeasurements | N | ?? |
| **OTHERS** ||||
| RepetitionTimes | RepetitionTimesExcitation | N | See extra 's', used(?) in case multiple TRs are passed -> not standard BIDS anyway |

Comments:
- **NAMING**. There are many other fields regarding the hardware, see the **Scanner hardware** section of [this file](https://github.com/bids-standard/bep001/blob/master/src/04-modality-specific-files/01-magnetic-resonance-imaging-data.md#scanner-hardware) and the same goes for the **Sequence Specifics**.
- **RF SPOILING**. There are 2 more fields *SpoilingState* (`true`/`false`) and *SpoilingType* (`RF`, `GRADIENT`, or `COMBINED`), see bottom of **Sequence Specifics** [section](https://github.com/bids-standard/bep001/blob/master/src/04-modality-specific-files/01-magnetic-resonance-imaging-data.md#sequence-specifics).
- **B1 MAPPING**. These are not part of the current proposal and would have to be added. What is missing is a clear definition of each of these parameters...
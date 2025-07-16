#ifndef INPUTFORM_H
#define INPUTFORM_H

#include <QWidget>
#include "OptionsData.h"

QT_BEGIN_NAMESPACE
class QLineEdit;
class QSpinBox;
class QDoubleSpinBox;
class QComboBox;
class QPushButton;
class QLabel;
QT_END_NAMESPACE

class InputForm : public QWidget
{
    Q_OBJECT

public:
    explicit InputForm(QWidget *parent = nullptr);

    // Get the current form data
    OptionsData getOptionsData() const;

    // Set form data
    void setOptionsData(const OptionsData &data);

    // Reset form to defaults
    void resetToDefaults();

signals:
    void advancedOptionsRequested();
    void runRequested();
    void formValidityChanged(bool isValid);

private slots:
    void onBrowseTPR();
    void onBrowseXTC();
    void onBrowseOutput();
    void onGridTypeChanged(int index);
    void validateForm();

private:
    void setupUI();
    void connectSignals();

    // Required fields
    QLineEdit *m_tprFileEdit;
    QLineEdit *m_xtcFileEdit;
    QPushButton *m_tprBrowseBtn;
    QPushButton *m_xtcBrowseBtn;

    // Grid size fields
    QComboBox *m_gridTypeCombo;
    QSpinBox *m_gridSingleSpin;
    QSpinBox *m_gridXSpin;
    QSpinBox *m_gridYSpin;
    QSpinBox *m_gridZSpin;
    QWidget *m_singleGridWidget;
    QWidget *m_tripleGridWidget;

    // Frame range
    QSpinBox *m_startFrameSpin;
    QSpinBox *m_endFrameSpin;
    QSpinBox *m_frameIntervalSpin;

    // Output
    QLineEdit *m_outputFileEdit;
    QPushButton *m_outputBrowseBtn;

    // Basic optional parameters
    QDoubleSpinBox *m_binSizeSpin;
    QDoubleSpinBox *m_qCutoffSpin;
    QComboBox *m_simulationTypeCombo;

    // Advanced options button
    QPushButton *m_advancedBtn;

    // Run button
    QPushButton *m_runBtn;

    // Status label
    QLabel *m_statusLabel;

    // Store advanced options
    OptionsData m_advancedOptions;
};

#endif // INPUTFORM_H
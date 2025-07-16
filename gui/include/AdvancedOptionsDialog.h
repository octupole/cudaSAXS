#ifndef ADVANCEDOPTIONSDIALOG_H
#define ADVANCEDOPTIONSDIALOG_H

#include <QDialog>
#include "OptionsData.h"

QT_BEGIN_NAMESPACE
class QSpinBox;
class QDoubleSpinBox;
class QLineEdit;
class QComboBox;
class QCheckBox;
QT_END_NAMESPACE

class AdvancedOptionsDialog : public QDialog
{
    Q_OBJECT

public:
    explicit AdvancedOptionsDialog(QWidget *parent = nullptr);

    void setOptionsData(const OptionsData &data);
    OptionsData getOptionsData() const;

private slots:
    void onScaledGridEnabledChanged(bool enabled);
    void onScaledGridTypeChanged(int index);
    void onWaterModelEnabledChanged(bool enabled);

private:
    void setupUI();
    void connectSignals();

    // B-spline order
    QSpinBox *m_orderSpin;

    // Scaled grid
    QCheckBox *m_scaledGridCheck;
    QComboBox *m_scaledGridTypeCombo;
    QSpinBox *m_scaledGridSingleSpin;
    QSpinBox *m_scaledGridXSpin;
    QSpinBox *m_scaledGridYSpin;
    QSpinBox *m_scaledGridZSpin;
    QDoubleSpinBox *m_scaleFactorSpin;
    QWidget *m_scaledGridWidget;
    QWidget *m_singleScaledWidget;
    QWidget *m_tripleScaledWidget;

    // Water model / padding
    QCheckBox *m_waterModelCheck;
    QWidget *m_waterModelWidget;
    QLineEdit *m_waterModelEdit;
    QSpinBox *m_sodiumSpin;
    QSpinBox *m_chlorineSpin;
};

#endif // ADVANCEDOPTIONSDIALOG_H